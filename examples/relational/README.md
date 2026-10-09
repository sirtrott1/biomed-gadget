# Cohort design analysis

`gadget_relational` analyzes a declared observational design complex and its
supported measurement fields. It retains source identities, explicit metrics,
boundary maps and complete fields through native homology and Hodge operations.

```sh
gadget_relational declaration.json build/cohort-results cohort_design_v1 1
```

The declaration must supply its actual carriers, incidence, measured fields,
support and source mappings. Equal dimensions do not identify carriers.
Missing observations require explicit support rather than replacement by zero.

The current profile uses positive unit patient and face metrics, relation
cardinality metrics, integral incidence and a Green parameter of one. These are
declared model choices. Observational comparisons do not establish causality or
predictive performance.

See the [C++ interface](../../include/gadget/relational.hpp) and
[scientific examples](../scientific/README.md) for the general application API.
