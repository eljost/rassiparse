<!doctype html>
<html>
    <head>
        <meta charset="utf-8">
        <style>
            html * {
                font-size: 1.1em;
            }
            .state {
                margin-bottom: 2%;
            }
            .mo {
                height: 10%;
                min-height: 5em;
                width: 10%;
                min-width: 5em;
                border: 1px solid #000000;
                padding: 1%;
            }
            sub {
                font-size: .8em;
            }
        </style>
    </head>
    <body>
    {% for sfs in sf_states %}
    <div class="state">
        <h2>
            {{ sfs.mult_label }}<sub>{{ sfs.state_rel }}</sub>,
            State {{ sfs.state }},
            Symmetry = {{ sfs.sym }},
            λ = {{ "%.1f" | format(sfs.dE_gs_nm) }} nm,
            ΔE = {{ "%.2f" | format(sfs.dE_gs_eV) }} eV,
            f = {{ "%.4f" | format(sfs.osc) }}
        </h2>
        {% for cdi in sfs.confdiff_images %}
        <div>
            <figure>
                <figcaption>{{ "%.1f" | format(cdi[0]*100) }}%</figcaption>
                {% for from_mo, to_mo in cdi[1] %}
                <img class="mo" src="{{ from_mo }}" />
                <svg width="7em" height="4em">
                  <defs>
                    <marker id="arrow" markerWidth="10" markerHeight="10" refX="0" refY="3" orient="auto" markerUnits="strokeWidth">
                      <path d="M0,0 L0,6 L9,3 z" fill="#000" />
                    </marker>
                  </defs>

                  <line x1="0" y1="25" x2="150" y2="25" stroke="#000" stroke-width="3" marker-end="url(#arrow)" />
                </svg>
                <img class="mo" src="{{ to_mo }}" />
                <br />
                {% endfor %}
            </figure>
        </div>
        {% endfor %}
        </div>
        {% endfor %}
    </body>
</html>
