#!/usr/bin/env python3

import argparse
import base64
import glob
import logging
import pprint
import os
import sys
import tarfile

import dash
import dash_core_components as dcc
import dash.dependencies as dep
import dash_html_components as html
from natsort import natsorted
import plotly.graph_objs as go

import numpy as np
import pandas as pd
import yaml

logging.basicConfig(level=logging.INFO)


def preload_pngs(png_fns):
    """Accepts a list of png filenames and returns a list
    of base64 encoded png images."""
    enc_pngs = list()
    for fn in png_fns:
        with open(fn, "rb") as handle:
            png = handle.read()
        enc_pngs.append(base64.b64encode(png))
        logging.info(f"Loaded '{fn}'.")

    return enc_pngs


def get_enc_src(enc_png):
    """Accepts a base64 encoded png and returns a string that can
    be used for the src attribute inside an img tag."""
    base_src = "data:image/png;base64,{}"
    return base_src.format(enc_png.decode())


def load_ressources(conf):
    energy_titles, energy_fns = zip(*[(title, conf["energies"][title])
                                      for title in conf["energies"]])
    energy_list = [(title, pd.read_csv(fn))
               for title, fn in zip(energy_titles, energy_fns)
               if os.path.exists(fn)
    ]

    # Load active space images
    active_space_glob = conf["active_space_glob"]
    as_fns = natsorted([fn for fn in glob.glob(active_space_glob)])
    active_space_imgs = preload_pngs(as_fns)

    # Load Mayer bond orders
    bo_glob = conf["bo_glob"]
    bo_df_fns = natsorted([fn for fn in glob.glob(bo_glob)])
    bo_dfs = [pd.read_csv(fn) for fn in bo_df_fns]
    for bo_df_fn in bo_df_fns:
        logging.info(f"Loaded '{bo_df_fn}'.")

    with open(conf["rassiscan"]) as handle:
        fns_and_weights = yaml.load(handle.read())

    # Preload individual MO images
    mo_glob = conf["mo_glob"]
    mo_fns = glob.glob(mo_glob)
    enc_mo_pngs = preload_pngs(mo_fns)
    mo_png_dict = {mo_fn: enc_mo_png
                   for mo_fn, enc_mo_png in zip(mo_fns, enc_mo_pngs)
    }

    return energy_list, active_space_imgs, bo_dfs, fns_and_weights, mo_png_dict


def pack(conf, conf_yaml):
    globbed = [conf_yaml, ]
    for to_glob in ("active_space_glob", "bo_glob", "mo_glob"):
        globbed.extend(glob.glob(conf[to_glob]))
    globbed.append(conf["rassiscan"])
    globbed.extend([conf["energies"][title] for title in conf["energies"]])

    packed_fn = "packed.tar.gz"
    with tarfile.open(packed_fn, "w:gz") as tar:
        for fn in globbed:
            tar.add(fn)
    logging.info(f"Packed files into {packed_fn}.")


def make_confdiff(enc_png_pairs, weight, mo_size):
    mo_width, mo_height = mo_size
    children = [html.H4("{:.1%}".format(weight))]
    for enc_from_png, enc_to_png in enc_png_pairs:
        from_src = get_enc_src(enc_from_png)
        to_src = get_enc_src(enc_to_png)
        children.extend([
                html.Img(
                    src=from_src,
                    width=mo_width,
                    height=mo_height
                ),

                html.Img(
                    src=to_src,
                    width=mo_width,
                    height=mo_height
                ),
            ]
        )
    return html.Div(children,  className="two columns")


def prepare_app(app, conf):
    ens_list, as_imgs, bo_dfs, fns_and_weights, mo_png_dict = load_ressources(conf)

    as_width, as_height = conf["active_space_size"]

    app.layout = html.Div([

        html.Div([
            # Potential energy curves
            html.Div([
                dcc.Dropdown(
                    id = "energies-selector",
                    options = [
                        {"label": title, "value": i}
                         for i, (title, ens_df) in enumerate(ens_list)
                    ],
                ),
                dcc.Graph(id = "energies",)],
                className = "six columns"
            ),

            # Active space images
            html.Div(
                html.Img(
                    id="active-space",
                    width=as_width,
                    height=as_height,
                ),
                className = "six columns"
            )],
            className = "row",
        ),

        html.Div([
            # Mayer Bond Orders
            html.Div(
                dcc.Graph(id="bond-orders"),
                className = "six columns"
            ),

            # Important configurations for a given step and root
            html.Div(id="fns-and-weights")
            ],
            className = "row",
        ),
    ])

    app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})


    @app.callback(
        dep.Output("energies", "figure"),
        [dep.Input("energies-selector", "value")]
    )
    def update_energies(value):
        if value is None:
            value = 0

        title, energies_df = ens_list[value]
        return {
            "data": [
                {"x": energies_df.index.values,
                 "y": energies_df[i],
                 "type": "lines",
                 "name": f"Root {int(i)+1}"
                } for i in energies_df
            ],
            "layout": go.Layout(
                title = title,
                xaxis = dict(
                    title = "Step"
                ),
                yaxis = dict(
                    title = "E / eV"
                ),
                hovermode = "closest"
            )
        }

    @app.callback(
        dep.Output("active-space", "src"),
        [dep.Input("energies", "hoverData")]
    )
    def update_active_space(hoverData):
        if hoverData is None:
            step_id = 0
        else:
            step_id = hoverData["points"][0]["x"]
        return "data:image/png;base64,{}".format(as_imgs[step_id].decode())

    @app.callback(
        dep.Output("bond-orders", "figure"),
        [dep.Input("energies", "clickData")]
    )
    def update_bond_orders(clickData):
        if clickData is None:
            root = 0
        else:
            root = clickData["points"][0]["curveNumber"]
        
        return {
            "data": [
                {"x": bo_dfs[root].index.values,
                 "y": bo_dfs[root][col],
                 "name": col
                } for col in bo_dfs[root]
            ],
            "layout": go.Layout(
                title = f"Mayer Bond Order - Root {root+1}",
                xaxis = dict(
                    title = "Step"
                ),
                yaxis = dict(
                    title = "Mayer BO"
                )
            )
        }


    @app.callback(
        dep.Output("fns-and-weights", "children"),
        [dep.Input("energies", "clickData")]
    )
    def update_fns_and_weights(clickData):
        if clickData is None:
            root = 0
            step = 0
        else:
            point = clickData["points"][0]
            step = point["pointNumber"]
            root = point["curveNumber"]

        cds = [
            html.H3(f"Step {step} -- Root {root+1}")
        ]
        for fns, weight in fns_and_weights[step][root]:
            enc_png_pairs = list()
            for from_fn, to_fn in fns:
                enc_from_png = mo_png_dict[from_fn]
                enc_to_png = mo_png_dict[to_fn]
                enc_png_pairs.append((enc_from_png, enc_to_png))

            cd = make_confdiff(enc_png_pairs, weight, conf["mo_size"])
            cds.append(cd)
            
        
        return cds


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("yaml")
    parser.add_argument("--port", default=8050, type=int)
    parser.add_argument("--pack", action="store_true")

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    with open(args.yaml) as handle:
        conf = yaml.load(handle.read())

    if args.pack:
        pack(conf, args.yaml)
        sys.exit()
    
    app = dash.Dash()
    prepare_app(app, conf)
    app.run_server(port=args.port)


if __name__ == "__main__":
    run()
