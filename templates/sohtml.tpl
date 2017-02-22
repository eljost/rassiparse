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
    {% for state in states %}
        <div class="state">
            <h2>
                {{ state.mult_label }}<sub>SO,{{ state.id }}</sub>,
                {# {{ state.sym }} #} Sym einbauen!,
                {{ "%.1f" | format(state.Enm) }} nm,
                {{ "%.2f" | format(state.EeV) }} eV,
                f={{ "%.4f" | format(state.osc) }}
            </h2>
            {% for from_img, to_img, weight in sf_trans_dict[state.sf_state] %}
            <div>
                <figure>
                    <img class="mo" src="{{ from_img }}" />
                    <img class="mo" src="{{ to_img }}" />
                    <figcaption>{{ "%.1f" | format(weight*100) }}%</figcaption>
                </figure>
            </div>
        {% endfor %}
        </div>
    {% endfor %}
    </body>
</html>
