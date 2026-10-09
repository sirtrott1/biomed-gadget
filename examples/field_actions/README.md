# Factored exact Green actions

Build the examples and run `gadget_field_actions_example 65536`. The argument
sets the number of retained coordinates and incidence occurrences. The authored
complex has one declared signed relation with coefficient column
`b = (1, -1/(n-1), ..., -1/(n-1))` and identity grade metrics. The source is the
first coordinate vector.

The example calls `NativeStudy::green_action` and checks every output coordinate
against the independent rank one identity
`(I + b b*)^-1 j = j - b(b*j)/(1 + b*b)`. It also verifies the complete native
operator residual. Source preparation, Green execution and verification timings
are separate. The result reports the two symmetric coordinate values and their
sum after computing and checking the complete field.

This operator has two eigenspaces and the exact action terminates in two steps.
Its retained boundary has n occurrences, while a materialized outer product can
have n squared entries. It demonstrates scaling through the representation and
native action contract. It does not establish the same solve time for an arbitrary
connected complex or compare mathematical analysis with neural model training.
