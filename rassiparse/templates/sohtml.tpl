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
                height: 15%;
                min-height: 10em;
                width: 15%;
                min-width: 10em;
                border: 1px solid #000000;
                padding: 1%;
            }
            sub {
                font-size: .8em;
            }
        </style>
    </head>
    <body>
    {% for so_state in so_states %}
        <div class="state">
            <h2>
                {{ so_state.mult }}<sub>SO{{ so_state.sostate }}</sub>,
                {{ "%.1f" | format(so_state.dE_global_nm) }} nm,
                {{ "%.2f" | format(so_state.dE_global_eV) }} eV,
                f={{ "%.5f" | format(so_state.osc) }}
            </h2>

        {% for sfs, weight in so_state.sf_states %}
            {% for cdi in sfs.confdiff_images %}
            <div>
                <figure>
                    <figcaption>{{ "%.1f" | format(cdi[0]*100) }}%<br />
                    {{ cdi[2].conf }}</figcaption>
                    {% for from_mo, to_mo in cdi[1] %}
                    <img class="mo" src="{{ from_mo }}" />
                    <svg width="7em" height="4em">
                      <defs>
                        <marker id="arrow" markerWidth="10" markerHeight="10" refX="0"
                                refY="3" orient="auto" markerUnits="strokeWidth">
                          <path d="M0,0 L0,6 L9,3 z" fill="#000" />
                        </marker>
                      </defs>

                      <line x1="0" y1="15" x2="110" y2="15" stroke="#000"
                            stroke-width="3" marker-end="url(#arrow)" />
                    </svg>
                    <img class="mo" src="{{ to_mo }}" />
                    <br />
                    {% endfor %}
                </figure>
                <br />
            </div>
            {% endfor %}
        {% endfor %}
        </div>
    {% endfor %}
    </body>
</html>
