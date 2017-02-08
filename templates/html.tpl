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
    {% for state in states[1:] %}
    <div class="state">
        <h2>
            S<sub>{{ state.id }}</sub>,
            {{ state.sym }},
            {{ "%.1f" | format(state.Enm) }} nm,
            {{ "%.2f" | format(state.EeV) }} eV,
            f={{ "%.4f" | format(state.f) }}
        </h2>
        {% if state.key in verbose_confs %}
        {% for start, end, weight in verbose_confs[state.key] %}
        <div>
            <figure>
                <img class="mo" src="{{ imgs[state.key[0]][start] }}" />
                <img class="mo" src="{{ imgs[state.key[0]][end] }}" />
                <figcaption>{{ "%.1f" | format(weight*100) }}%</figcaption>
            </figure>
        </div>
        {% endfor %}
        {% endif %}
        </div>
        {% endfor %}
    </body>
</html>
