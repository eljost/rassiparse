{% extends "rassi_base.tpl" %}
{% block visualization %}
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

          <line x1="0" y1="15" x2="110" y2="15" stroke="#000" stroke-width="3" marker-end="url(#arrow)" />
        </svg>
        <img class="mo" src="{{ to_mo }}" />
        <br />
        {% endfor %}
    </figure>
</div>
{% endfor %}
{% endblock %}
