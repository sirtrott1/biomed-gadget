# Prepared temporal queries

A `ScientificSession` reuses prepared operators across queries on the same
admitted study. Each query still declares its source endpoint, fields and
metrics; preparations are reused only when those inputs agree.

```cpp
gadget::ScientificSession session(study);
auto first = session.execute_record(first_query);
auto second = session.execute_record(second_query);
```

The returned records retain the requests and native result fields. Keep their
owners alive when accessing field or tensor references. The
[scientific sequence example](../../scientific/README.md) runs an ordered batch
through this interface, and the [C++ guide](../../../include/gadget/README.md)
describes ownership and direct library use.
