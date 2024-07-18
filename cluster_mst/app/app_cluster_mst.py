# panel serve <filename>
import base64
from io import StringIO
from io import BytesIO

import pandas as pd
# import numpy as np

# import altair as alt
import holoviews as hv
from holoviews.streams import Selection1D
hv.extension("bokeh")
from bokeh.models import HoverTool

import panel as pn
import panel.widgets as pnw

from jupy_tools import utils as u
from jupy_tools import mol_view as mv, cluster_mst as cmst

pn.extension()

HELP_TEXT = """
# Cluster MST

This app takes a file containing structures as Smiles and activity data, and generates a clustering of the structures, represented as Minimum Spanning Tree (MST).
The input file has to be a &lt;tab&gt;-separated TSV file and contain at least three columns: &lt;Identifier&gt; (default: Compound_Id), &lt;Activity&gt; (check the "Reverse" box if lower values are better, e.g. for "Activity"), and "Smiles" for the structures.
The tool takes the top N active compounds ("Top N active", default: 50) and adds the most similar compounds for each of these ("Number of similar compounds", default: 10), downto a minimum similarity cutoff ("Similarity cutoff", default: 0.6) using the chosen fingerprint method, then generates the MST.

Generally, only linear-scaled values should be used for Activity, e.g. percentages. IC50 values should be converted to pIC50.
The points in the plot can be selected using the lasso tool, and the selected compounds are shown in the table below the plot.
"""

button = pnw.Button(name="Start Calculation", button_type="primary", disabled=False)

def mol_image_tag(mol):
    mol_tag = mv.MolImage(mol, svg=False, size=250).tag  # options='width="70%"'
    return mol_tag


def struct_hover(id_col, cols):
    """Create a structure tooltip that can be used in Holoviews.
    Takes a MolFrame instance as parameter."""
    add_cols = []
    if cols is None:
        cols = []
    if isinstance(cols, str):
        cols = [cols]
    for col in cols:
        add_cols.append(
            f"""<div>
                    <span style="font-size: 12px;">{col}: @{col}</span>
                </div>"""
        )
    add_cols_txt = "\n".join(add_cols)
    # <img src="@Image" alt="Mol" width="70%"><br>
    hover = HoverTool(
        tooltips=f"""
            <div>
                <div>
                    @Image<br>
                <div>
                <div>
                    <span style="font-size: 12px; font-weight: bold;">@{id_col}</span>
                </div>
                {add_cols_txt}
            </div>
        """
    )
    return hover


@pn.depends(button)
def show_result(event=None):    
    def update_table(index):
        selected_df = mst.df.iloc[index].copy()
        result[1] = pn.pane.DataFrame(selected_df[columns], escape=False, index=False, max_width=500)
        sio = StringIO()
        selected_df[dl_columns].to_csv(sio, sep="\t", index=False)
        sio.seek(0)
        result[2] = pnw.FileDownload(sio, embed=True, filename='selection.tsv')


    print("In function")
    print(w_file_input.filename)
    if w_file_input.filename is None:
        return pn.pane.Markdown(HELP_TEXT + "\n\nPlease upload a file.")
    df = pd.read_csv(BytesIO(w_file_input.value), sep="\t")
    print(len(df))
    avail_cols = ", ".join(df.columns.tolist())
    if w_id_col.value not in df.columns:
        return pn.pane.Markdown(
            (
                HELP_TEXT + f"\n\nERROR: Identifier column {w_id_col.value} not found in the file."
                + f"\n\nAvailable columns: \n{avail_cols}"
            )
        )
    if w_act_col.value not in df.columns:
        return pn.pane.Markdown(
            (
                HELP_TEXT + f"\n\nERROR: Activity column {w_act_col.value} not found in the file."
                + f"\n\nAvailable columns: \n{avail_cols}"
            )
        )
    if "Smiles" not in df.columns:
        return pn.pane.Markdown(
            (
                HELP_TEXT + "\n\nERROR: Smiles column not found in the file."
                + f"\n\nAvailable columns: \n{avail_cols}"
            )
        )
    if w_top_n_act.value < 1:
        return pn.pane.Markdown(
            HELP_TEXT + "\n\nERROR: Top N active must be at least 1."
        )
    if w_top_n_act.value > len(df):
        w_top_n_act.value = len(df)
    if w_num_sim.value < 0:
        return pn.pane.Markdown(
            HELP_TEXT + "\n\nERROR: Number of similar compounds must be at least 0."
        )
    if w_sim_cutoff.value < 0.2 or w_sim_cutoff.value > 0.9:
        return pn.pane.Markdown(
            HELP_TEXT + "\n\nERROR: Similarity cutoff must be between 0.2 and 0.9."
        )
    
    dl_columns = df.columns.tolist()    
    mst = cmst.ClusterMST(
        df, id_col=w_id_col.value, act_col=w_act_col.value, 
        top_n_act=w_top_n_act.value, num_sim=w_num_sim.value, 
        reverse=w_reverse.value, sim_cutoff=w_sim_cutoff.value,
        fp=w_fp_method.value
    ) 
    mst.calc_mst()
    mst.df = u.calc_from_smiles(mst.df, "Image", mol_image_tag, smiles_col="Smiles")
    columns = ["Image", mst.id_col, mst.act_col]

    tooltip = [mst.act_col]
    edges = [[(x1, y1), (x2, y2)] for (x1, y1, x2, y2) in mst.edges.to_records(index=False)]

    hover = struct_hover(mst.id_col, cols=tooltip)
    plot_options = {
        "width": 1200,
        "height": 800,
        "tools": [hover, "lasso_select"],
        "legend_position": "right",
        "toolbar": "right",
        "size": 12,
        "colorbar": True,
    }
    colorby = mst.act_col
    kdims = ["X", "Y"]
    vdims = [mst.id_col, "Image", colorby]
    df = df.sort_values(colorby, ascending=mst.reverse)
    scatter = hv.Points(data=mst.df, kdims=kdims, vdims=vdims)  # , label=title)
    plot_options["color"] = colorby
    plot_options["cmap"] = w_cmap.value
    
    selection = Selection1D(source=scatter)
    selection.param.watch(lambda event: update_table(event.new), 'index')
    
    # chart = (hv.Path(edges).options(color="black") * scatter.options(**plot_options))
    chart = hv.Path(edges).opts(color="black") * scatter.opts(**plot_options)

    result = pn.Column(
        chart,
        pn.pane.Markdown("Select points in the chart using the lasso select tool."),
        pn.pane.Markdown(""),
    )

    return result


title = "Cluster MST"
w_file_input = pnw.FileInput(
    accept='.tsv', multiple=False,
    name="Upload a file", description="Select a file with structures and activity data."
)
w_id_col = pnw.TextInput(
    name="Identifier column", value="Compound_Id",
    description="Name of the column containing the compound identifiers."
)
w_act_col = pnw.TextInput(
    name="Activity column",
    description="Name of the column containing the activity values."
)
w_reverse = pnw.Checkbox(
    name="Reverse",
    # description="Check if lower values are better, e.g. for IC50."
)
w_top_n_act = pnw.IntInput(
    name="Top N active", value=50,
    description="Number of top active compounds to select."
)
w_num_sim = pnw.IntInput(
    name="Number of similar compounds", value=10,
    description="Number of similar compounds to select for each active compound."
)
w_sim_cutoff = pnw.FloatInput(
    name="Similarity cutoff", value=0.6,
    start=0.2, end=0.9, step=0.1,
    description="Minimum similarity cutoff."
)
w_fp_method = pnw.Select(
    name="Fingerprint method",
    options=sorted(cmst.FPDICT.keys()),
    value="ECFC4",
    description="Fingerprint method to use for similarity calculation."
)
w_cmap = pnw.Select(
    name="Color map",
    options=["brg", "bmy", "viridis", "plasma", "magma", "turbo"],
    value="brg",
    description="Color map for the plot."
)


print("Starting...")
app = pn.template.FastListTemplate(
    site="Datavis",
    title=title,
    theme_toggle=False,
    sidebar=[
        pn.pane.Markdown("### Upload file with structures and activity data."),
        w_file_input,
        w_id_col,
        w_act_col,
        w_reverse,
        w_top_n_act,
        w_num_sim, 
        w_sim_cutoff,
        w_fp_method,
        w_cmap,
        
        button
    ],
    main=[
        pn.Column(
            show_result,
        )
    ],
    main_layout=None,  # Removes the grey background around the results in the main panel
)

app.servable()
