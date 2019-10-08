{{ fullname | escape | underline }}

.. rubric:: Description

.. automodule:: {{ fullname }}

.. currentmodule:: {{ fullname }}


{% if classes %}
.. rubric:: Classes

.. autosummary::
    :toctree: .
    :nosignatures:
    :template: class.rst

    {% for class in classes %}
    {{ class }}
    {% endfor %}

{% endif %}


{% if functions %}
.. rubric:: Functions

.. autosummary::
    :toctree: .
    :nosignatures:

    {% for function in functions %}
    {{ function }}
    {% endfor %}

{% endif %}


