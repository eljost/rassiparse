#!/usr/bin/env python3

import argparse
import glob
import logging
import os
import pathlib
import sys

from jinja2 import Environment, PackageLoader
from natsort import natsorted
import yaml

logging.basicConfig(level=logging.INFO)


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("yaml")

    return parser.parse_args(args)

    parser.add_argument("steps", type=int)
    parser.add_argument("basis", default="ano-rcc-vdzp")
    parser.add_argument("shift", default="0.3")


def run():
    args = parse_args(sys.argv[1:])

    with open(args.yaml) as handle:
        conf = yaml.load(handle.read())

    # Load templates
    env = Environment(loader=PackageLoader('rassiparse', 'templates'))
    inp_tpl = env.get_template(conf["template"])
    sub_tpl = env.get_template("submolcas.tpl")
    run_tpl = env.get_template("runscan.tpl")

    coord_fns = natsorted(glob.glob(conf["coord_glob"]))
    jobiph_fns = natsorted(glob.glob(conf["jobiph_glob"]))

    assert(len(coord_fns) == len(jobiph_fns))

    basis = conf["basis"]
    states = " ".join([str(i) for i in range(1, conf["multistate"]+1)])
    multistate = "{} {}".format(conf["multistate"], states)
    shift = conf["shift"]
    root = pathlib.Path(conf["root_dir"])
    step_fmt = conf["step_fmt"]
    project = conf["project"]
    input_fn_base = conf["input_fn_base"]

    run_fn = root / "run.sh"

    step_strs = list()
    for step, (coord, jobiph) in enumerate(zip(coord_fns, jobiph_fns)):
        logging.info(f"Step {step}:")
        logging.info(f"{coord}, {jobiph}.")
        step_formatted = step_fmt.format(step)
        step_str = f"step{step_formatted}"
        step_strs.append(step_str)
        job_name = step_str

        # Generate curr_dir for this job
        curr_dir = root / f"step{step_formatted}"
        input_fn_no_ext = f"{input_fn_base}.step{step_formatted}"
        input_fn = f"{input_fn_no_ext}.in"
        input_fn_path = curr_dir / input_fn
        sub_fn_path = curr_dir / "submolcas.sh"

        # Render MOLCAS input
        inp_rendered = inp_tpl.render(
                            coord_fn=coord,
                            basis=basis,
                            jobiph_fn=jobiph,
                            multistate=multistate,
                            shift=shift,
                            curr_dir=curr_dir,
                            step_str=step_str,
        )

        # Render Submit script
        sub_rendered = sub_tpl.render(
                            job_name=job_name,
                            curr_dir=curr_dir,
                            project=project,
                            input_fn_no_ext=input_fn_no_ext,
                            step=step_formatted,
        )
        try:
            os.mkdir(curr_dir)
        except FileExistsError:
            pass

        with open(input_fn_path, "w") as handle:
            handle.write(inp_rendered)

        with open(sub_fn_path, "w") as handle:
            handle.write(sub_rendered)

    run_rendered = run_tpl.render(step_strs=step_strs)

    with open(run_fn, "w") as handle:
        handle.write(run_rendered)

if __name__ == "__main__":
    run()
