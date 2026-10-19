# Reusing preparations for temporal relational fields

Reuse prepared operators for an ordered sequence of field queries. Each result retains the source version and the maps needed to interpret its coordinates.

A scientific session reuses preparations across queries on the same admitted
source and declared metrics:

```cpp
gadget::ScientificSession session(study);
auto first = session.execute(query);
auto second = session.execute(query);
```

Use the [scientific query schema](../../../schemas/scientific-query.schema.json)
for the declared operations and the [C++ interface](../../../include/gadget/README.md)
for direct library calls. The [first analysis](../../../README.md#first-analysis)
shows a complete synthetic query through the command line application.
