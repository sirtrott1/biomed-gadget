"""Installed contract inspection keeps execution and qualification separate."""
import sys
import unittest
from unittest.mock import patch

from gadget_data.capabilities import native_primitive_capabilities, operation_capabilities


class Capabilities(unittest.TestCase):
    def test_unknown_operation_is_not_an_execution_contract(self):
        row, = operation_capabilities(("NO_SUCH_GADGET_OPERATOR",))
        self.assertEqual(row["status"], "unavailable")
        self.assertIsNone(row["contract"])

    def test_file_actions_preserve_permissions_and_effects(self):
        run, replay = operation_capabilities(("DEMAND_RUN", "DEMAND_REPLAY"))
        for row in (run, replay):
            if row["contract"] is None:
                self.assertEqual(row["status"], "unavailable")
                continue
            contract = row["contract"]["current"]
            self.assertTrue({"read", "files", "file_read"} <= set(contract["requires"]))
            self.assertFalse(contract["memoizable"])
            self.assertEqual(row["platform_qualification"], "unverified")
            self.assertEqual(row["biomedical_qualification"], "unverified")

    def test_constructor_observation_and_action_scopes_remain_separate(self):
        report = native_primitive_capabilities()
        if not report["available"]:
            self.assertEqual(report["primitives"], ())
            return
        rows = {row["name"]: row for row in report["primitives"]}
        self.assertEqual(rows["UnivariateFamily"]["inspection"], "declaration-only")
        self.assertEqual(rows["Domain"]["inspection"], "native-observation")
        self.assertEqual(rows["PreparedNativeMoments"]["status"], "native-unregistered")
        self.assertEqual(rows["ActionMoments"]["status"], "native-qualified")

    def test_numerical_and_native_contracts_are_not_aliased(self):
        rows = operation_capabilities(("GREEN", "GREEN_FIELD", "C_EVALUATE"))
        self.assertEqual([row["name"] for row in rows], ["GREEN", "GREEN_FIELD", "C_EVALUATE"])
        contracts = [row["contract"]["current"]["implementation_key"] for row in rows if row["contract"] is not None]
        self.assertEqual(len(contracts), len(set(contracts)))

    def test_inspection_never_calls_an_operator_or_loads_optional_services(self):
        before = {key for key in sys.modules if key.startswith(("agent", "rexgraph_cpp", "system"))}
        with patch("rcql.operators.get_operator", side_effect=AssertionError("execution lookup")):
            operation_capabilities(("DEMAND_RUN", "C_WORKFLOW_BUILD"))
        after = {key for key in sys.modules if key.startswith(("agent", "rexgraph_cpp", "system"))}
        self.assertEqual(before, after)

    def test_invalid_names_are_refused(self):
        for names in ((None,), (1,), ("",)):
            with self.assertRaises(TypeError):
                operation_capabilities(names)


if __name__ == "__main__":
    unittest.main()
