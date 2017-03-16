{% extends "rassi_base.tpl" %}
{% block visualization %}
{% for weight, mos, spins in sfs.single_mos %}
<div>
    <table>
        <tr>
            <th></th>
        {% for spin in spins %}<th>{{ spin }}</th>{% endfor %}
        </tr>
        <tr>
            <td>{{ "%.1f" | format(weight*100) }}%</td>
        {% for mo in mos %}<td><img class="singlemo" src="{{ mo }}" /></td>{% endfor %}
        </tr>
    </table>
</div>
{% endfor %}
{% endblock %}
