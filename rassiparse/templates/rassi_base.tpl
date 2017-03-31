<!doctype html>
<html>
    <head>
        <meta charset="utf-8">
        <style>
            html * {
                font-size: 1.05em;
            }
            .state {
                margin-bottom: 5%;
            }
            .mo {
                min-height: 5em;
                max-height: 7em;
                min-width: 7em;
                max-width: 9em;
                border: 1px solid #000000;
            }
            .singlemo {
                min-height: 3em;
                max-height: 5em;
                min-width: 5em;
                max-width: 8em;
                border: 1px solid #000000;
            }
            sub {
                font-size: .8em;
            }
            table {
                margin: 0;
            }
            h3 {
                margin: 0;
            }
            .weight_conf {
                margin-left: 2%;
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
    {% block visualization scoped %}{% endblock %}
    </div>
    {% endfor %}
    </body>
</html>
