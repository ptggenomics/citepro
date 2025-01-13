from ipyfilechooser import FileChooser
import ipywidgets as widgets
from IPython.display import display
import celltypist as ct

from anndata import AnnData

from .pl import multiomics_feature_plot


modname = [mn for mn in ct.models.models_description()['model']]
modname.insert(0, "--- None ---")



def generate_create_mu_setting():
    
    # Create and display a FileChooser widget
    fc_count = FileChooser()
    fc_count.filter_pattern = "*.h5"
    fc_count.show_only_dirs = True
    fc_count_reset = widgets.Button(description='Reset selection', layout=widgets.Layout(width='150px'))
    fc_count_reset.on_click(lambda x: fc_count.reset())

    #fc_count = 'Count_matrix'
    create_mu_setting = {}

    ta_samp_id = widgets.Text(
        description='Sample_ID'
    )
    # Create another FileChooser widget for barcode block list
    fc_barcode = FileChooser()
    fc_barcode.show_only_dirs = True
    fc_barcode_reset = widgets.Button(description='Reset selection', layout=widgets.Layout(width='150px'))
    fc_barcode_reset.on_click(lambda x: fc_barcode.reset())
    #fc_barcode.title = 'Barcode_block_list'

    prot_norm_dropdown = widgets.Dropdown(
        options=['asinh', 'clr', 'none'],
        value='asinh'
    )

    celltypist_model_dropdown = widgets.Dropdown(
        options=modname,
        value='Immune_All_Low.pkl'
    )

    cb_3d_umap = widgets.Checkbox(value=True, description='Generate 3D umap')

    def on_set(b):
        print('Check parameters...')
        if not fc_count.selected:
            print('Please select a count matrix.')
            return
        ct_model = celltypist_model_dropdown.value if celltypist_model_dropdown.value != "--- None ---" else None
        sampid = ta_samp_id.value if ta_samp_id.value else None
        
        print("Setting:\n    count_matrix: {} \n    barcode_allow_list: {}\n    sample_id: {}\n    prot_norm: {}\n    celltypist_model: {}\n    add_3d_umap: {}".format(
            fc_count.selected, fc_barcode.selected, sampid, prot_norm_dropdown.value, ct_model, cb_3d_umap.value))

        nonlocal create_mu_setting

        create_mu_setting['path_count'] = fc_count.selected
        create_mu_setting['allow_file'] = fc_barcode.selected
        create_mu_setting['samp_id'] = sampid
        create_mu_setting['prot_norm'] = prot_norm_dropdown.value
        create_mu_setting['celltypist_model'] = ct_model
        create_mu_setting['add_3d_umap'] = cb_3d_umap.value

        print('Good to go!')

    set_button = widgets.Button(description='Set parameters', layout=widgets.Layout(width='150px'))
    set_button.on_click(on_set)

    lay1 = widgets.Layout(border = '1px solid #777777', padding = '2px', margin = '5px')
    hboxes = []
    hboxes.append(widgets.HBox([widgets.Label(value = 'Count_matrix'), fc_count, fc_count_reset], layout = lay1))
    hboxes.append(widgets.HBox([widgets.Label(value = 'Barcode_allow_list'), fc_barcode, fc_barcode_reset], layout = lay1))
    hboxes.append(widgets.HBox([ta_samp_id], layout = lay1))
    hboxes.append(widgets.HBox([widgets.Label(value = 'Protein normalization method'), prot_norm_dropdown], layout = lay1))
    hboxes.append(widgets.HBox([widgets.Label(value = 'Celltypist model'), celltypist_model_dropdown], layout = lay1))
    hboxes.append(widgets.HBox([cb_3d_umap], layout = lay1))
    hboxes.append(widgets.HBox([set_button], layout = lay1))
    #hboxes.append(widgets.HBox(output, layout = lay1))

    #hbox2 = widgets.HBox([ta_samp_id, celltypist_model_dropdown])                     
    # Wrap the file choosers with VBox
    vbox = widgets.VBox(hboxes)
    display(vbox)

    return create_mu_setting

def display_multiomics_feature_plot(adata: AnnData):

    prot_idx = adata.var[adata.var['feature_types'] == 'Antibody Capture' ].index
    
    cb_prot = widgets.Combobox(
        # value='John',
        placeholder='prot:CD8A.65146.1',
        options=[pt for pt in prot_idx],
        description='Protein:',
        ensure_option=True,
        disabled=False
    )

    output = widgets.Output()

    def on_value_change(change):
        with output:
            output.clear_output()
            multiomics_feature_plot(adata, change.new)

    cb_prot.observe(on_value_change, names='value')

    display(cb_prot, output)