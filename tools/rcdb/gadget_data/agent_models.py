"""Hive recipe proposals admitted against GADGET's actual native declarations."""
from collections.abc import Mapping
from dataclasses import replace
import json


class GadgetRecipeWorker:
    """A Hive capability over pinned workspace/intent records and one model bee.

    The model chooses declared goal names or asks for clarification. It cannot
    supply numerical parameters, source bindings or an execution policy. Native
    GADGET admits the proposal; automatic execution requires the caller's typed
    answer contract. Complete admission families and result publications remain
    in the returned RCQL Result, including untyped and unsupported proposals.
    """

    def __init__(self, service, hive, model, *, workspace, intent, execute=True, generation=None):
        def pin(value):
            if (not isinstance(value, Mapping) or set(value) != {"record_id", "version"}
                    or type(value["record_id"]) is not str or not value["record_id"]
                    or type(value["version"]) is not int or value["version"] < 1):
                raise ValueError("recipe sources require explicit record_id/version pins")
            return value["record_id"], value["version"]
        if type(execute) is not bool:
            raise TypeError("recipe execute must be bool")
        options = dict(generation or {})
        if set(options) - {"max_tokens", "temperature", "chat_template_kwargs"}:
            raise ValueError("recipe generation accepts max_tokens, temperature and chat_template_kwargs")
        self.service, self.hive, self.model = service, hive, model
        self._workspace, self._intent = pin(workspace), pin(intent)
        self._execute = execute
        # Retain caller configuration independently of later dictionary mutation.
        self._generation = json.dumps(options, allow_nan=False)

    def __call__(self, arguments):
        from rcql import call, query, source
        if (not isinstance(arguments, Mapping) or "question" not in arguments
                or set(arguments) - {"question", "save_as", "expected_version"}
                or type(arguments["question"]) is not str or not arguments["question"]):
            raise ValueError("recipe worker requires question and optional checked result publication")
        question = arguments["question"]
        pins = (*self._workspace, *self._intent)
        catalog_result = self.service.query(query(source("db"), call("GADGET_RECIPE_CATALOG", *pins)))
        catalog = json.loads(catalog_result.values[0])
        names = [entry["id"] for entry in catalog["recipes"]]
        shared = {"schema": {"const": "gadget.pipeline.recipe.proposal.v1"},
                  "source_digest": {"const": catalog["source_digest"]},
                  "intent_sha256": {"const": catalog["intent_sha256"]}}
        clarify = {"type": "object", "properties": {
            "action": {"const": "clarify"}, "reason": {"type": "string", "minLength": 1}, **shared},
            "required": ["action", "reason", *shared], "additionalProperties": False}
        branches = [clarify]
        if names:
            branches.append({"type": "object", "properties": {
                "action": {"const": "recipes"}, "goals": {"type": "array", "minItems": 1,
                    "maxItems": len(names),
                    "items": {"type": "string", "enum": names}}, **shared},
                "required": ["action", "goals", *shared], "additionalProperties": False})
        prompt = json.dumps({"task": "Choose declared goal bindings that answer this request. "
            "Use clarify when the request is unsupported, ambiguous, negated or lacks required inputs. "
            "The catalog contains actual supplied declarations; do not invent parameters or data. "
            "A named operation is not evidence that it answers the question. "
            "Return a GADGET recipe proposal using the declared JSON schema.",
            "catalog": catalog, "request": question}, ensure_ascii=False)
        reply = self.hive.ask_full(self.model, prompt,
            system="Propose GADGET recipes bound to the supplied sources or request clarification. Native GADGET checks every proposal.",
            response_format={"type": "json_schema", "json_schema": {"name": "gadget_recipe_proposal",
                "strict": True, "schema": {"oneOf": branches}}}, **json.loads(self._generation))
        proposal = None
        unresolved = "model_unavailable" if reply is None else (
            "incomplete_model_output" if reply.finish_reason == "length" else None)
        if unresolved is None:
            try:
                proposal = json.loads(reply.content or "")
            except json.JSONDecodeError:
                unresolved = "invalid_model_output"
        # Native admission validates the complete schema, names, actual sources,
        # original question and exact finite answer obligations before execution.
        if unresolved is None:
            result = self.service.query(query(source("db"), call("GADGET_RECIPE_PROPOSAL", *pins,
                reply.content, question, self._execute)))
        else:
            from rcql import Exactness, Result
            result = Result(({"observation": {"status": unresolved, "executed": False},
                              "publications": {}},), plan=("GADGET_RECIPE_UNRESOLVED",),
                exactness=(Exactness.STRUCTURAL,), provenance=catalog_result.provenance,
                execution=catalog_result.execution)
        result = replace(result, native_plan={**(result.native_plan or {}),
            "gadget_recipe_model": {"model": self.model, "question": question,
                "proposal": proposal, "raw_reply": None if reply is None else reply.raw,
                "catalog_execution": catalog_result.execution,
                "catalog_packet": catalog_result.to_bytes()}})
        if arguments.get("save_as") is not None:
            self.service.session.save_result(arguments["save_as"], result,
                expected_version=arguments.get("expected_version", 0))
        return result
