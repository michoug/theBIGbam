#!/usr/bin/env python3
"""
Start a small Bokeh server that shows a collapsible control panel and a main panel.
The control panel is populated from the SQLite database (`Sample`, `Contig`, `Variable` tables).
When the Apply button is clicked the server calls `generate_bokeh_plot` from `plotting_data.py`
and replaces the main panel with the returned Bokeh layout.

Run:
    python start_bokeh_server.py --db path/to/db --port 5006

Or just run and open the printed URL.
"""

import argparse
import sqlite3
import traceback

from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import Div
from bokeh.models.widgets import Select, CheckboxGroup, CheckboxButtonGroup, Button

# Import the plotting function from the repo
from plotting_data import generate_bokeh_plot


def build_controls(conn):
    """Query DB and return widgets and helper mappings."""
    cur = conn.cursor()

    # Samples
    cur.execute("SELECT Sample_name FROM Sample ORDER BY Sample_name")
    samples = [r[0] for r in cur.fetchall()]
    if not samples:
        samples = [""]
    sample_select = Select(title="Sample", value=samples[0], options=samples)

    # Contigs
    cur.execute("SELECT Contig_name FROM Contig ORDER BY Contig_name")
    contigs = [r[0] for r in cur.fetchall()]
    if not contigs:
        contigs = [""]
    contig_select = Select(title="Contig", value=contigs[0], options=contigs)

    # Modules and variables
    cur.execute("SELECT DISTINCT Module FROM Variable ORDER BY Module")
    modules = [r[0] for r in cur.fetchall()]
    if not modules:
        modules = []

    module_checkbox = CheckboxGroup(labels=modules, active=list(range(len(modules))))

    # For each module get variables
    module_to_vars = {}
    module_var_widgets = {}
    for i, mod in enumerate(modules):
        cur.execute("SELECT Variable_name FROM Variable WHERE Module=? ORDER BY Variable_name", (mod,))
        vars_for_mod = [r[0] for r in cur.fetchall()]
        module_to_vars[mod] = vars_for_mod
        if len(vars_for_mod) > 1:
            # use CheckboxGroup for selecting individual variables
            w = CheckboxGroup(labels=vars_for_mod, active=list(range(len(vars_for_mod))))
            module_var_widgets[mod] = w
        else:
            # no widget needed for single-variable modules
            module_var_widgets[mod] = None

    apply_button = Button(label="Apply", button_type="primary")

    widgets = {
        'sample_select': sample_select,
        'contig_select': contig_select,
        'module_checkbox': module_checkbox,
        'module_var_widgets': module_var_widgets,
        'apply_button': apply_button,
        'modules': modules,
        'module_to_vars': module_to_vars,
    }
    return widgets


def modify_doc_factory(db_path, max_visible_width=1800, subplot_size=130):
    """Return a modify_doc(doc) function to be used by Bokeh server application."""
    conn = sqlite3.connect(db_path)

    widgets = build_controls(conn)

    instructions = Div(text="<b>Controls</b>: select sample, contig and modules/variables then click Apply.")

    controls_children = [instructions, widgets['sample_select'], widgets['contig_select'], widgets['module_checkbox']]
    # append variable selectors
    for mod in widgets['modules']:
        var_widget = widgets['module_var_widgets'].get(mod)
        if var_widget is not None:
            # label
            controls_children.append(Div(text=f"<b>{mod} variables</b>"))
            controls_children.append(var_widget)

    controls_children.append(widgets['apply_button'])

    controls_column = column(*controls_children, width=350)
    # Try to use Panel/Accordion when available (Bokeh versions vary). Otherwise fall back to plain column.
    try:
        from bokeh.models.widgets import Panel, Accordion  # import here to avoid import-time error
        accordion_panel = Panel(child=controls_column, title="Controls")
        accordion = Accordion(children=[accordion_panel])
    except Exception:
        accordion = controls_column

    main_placeholder = column(Div(text="<i>No plot yet. Select options and click Apply.</i>"))

    layout = row(accordion, main_placeholder)

    def apply_clicked():
        try:
            sample = widgets['sample_select'].value
            contig = widgets['contig_select'].value

            # Build requested_features list
            requested = []
            modules = widgets['modules']
            active_modules = set(widgets['module_checkbox'].active)
            for idx, mod in enumerate(modules):
                vars_for_mod = widgets['module_to_vars'].get(mod, [])
                var_widget = widgets['module_var_widgets'].get(mod)
                if idx in active_modules:
                    # module selected -> include all module variables
                    requested.extend(vars_for_mod)
                else:
                    # module not selected -> if var widget exists, include individually selected vars
                    if var_widget is not None:
                        for vidx in var_widget.active:
                            requested.append(var_widget.labels[vidx])

            # Deduplicate while preserving order
            seen = set()
            requested_features = [x for x in requested if not (x in seen or seen.add(x))]
            if not requested_features:
                main_placeholder.children = [Div(text="<b>No variables selected.</b>")]
                return

            # Call plotting_data.generate_bokeh_plot to build the grid
            grid = generate_bokeh_plot(db_path, requested_features, contig, sample, max_visible_width, subplot_size)

            # Replace main panel
            main_placeholder.children = [grid]

        except Exception as e:
            tb = traceback.format_exc()
            main_placeholder.children = [Div(text=f"<pre>Error building plot:\n{tb}</pre>")]

    widgets['apply_button'].on_click(lambda: apply_clicked())

    def modify_doc(doc):
        doc.add_root(layout)

    return modify_doc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True, help="Path to sqlite DB")
    parser.add_argument("--port", type=int, default=5006, help="Port to serve Bokeh app")
    parser.add_argument("--width", type=int, default=1800, help="Plot width")
    parser.add_argument("--height", type=int, default=130, help="Subplot height")
    args = parser.parse_args()

    modify_doc = modify_doc_factory(args.db, max_visible_width=args.width, subplot_size=args.height)

    # Start bokeh server programmatically
    try:
        from bokeh.server.server import Server
        from bokeh.application import Application
        from bokeh.application.handlers.function import FunctionHandler
    except Exception as e:
        print("Bokeh server components not installed or unavailable:", e)
        print("You can run this file with `bokeh serve --show start_bokeh_server.py --args --db <DB>` as alternative.")
        return

    app = Application(FunctionHandler(modify_doc))
    server = Server({'/': app}, port=args.port)
    server.start()
    print(f"Bokeh server running at http://localhost:{args.port}/")
    server.io_loop.add_callback(server.show, "/")
    server.io_loop.start()


if __name__ == "__main__":
    main()
