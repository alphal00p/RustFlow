#!/usr/bin/env python3
"""Focused launcher contract tests; no Wolfram installation or license needed."""

import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


source = Path(__file__).with_name("wolfram_one_thread.py")
spec = importlib.util.spec_from_file_location("wolfram_one_thread", source)
launcher = importlib.util.module_from_spec(spec)
spec.loader.exec_module(launcher)


class LauncherTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="wolfram launcher ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.kernel = self.directory / "licensed kernel"
        self.kernel.write_text("#!/bin/sh\nexit 0\n")
        self.kernel.chmod(0o755)
        self.script = self.directory / 'original "script".wl'
        self.script.write_text("Exit[0];\n")
        self.environment = patch.dict(
            os.environ, {"RUSTFLOW_WOLFRAM_KERNEL": str(self.kernel)}, clear=True
        )
        self.environment.start()
        self.addCleanup(self.environment.stop)

    def test_flags_script_arguments_and_exact_filenames_are_preserved(self):
        arguments = ["-noinit", "-noprompt", "-script", str(self.script),
                     "with spaces", 'a"b', "literal$(command)", "-script", "λ"]
        with patch.object(launcher.os, "execv") as execute:
            launcher.launch(arguments)
        expected = arguments.copy()
        expected[3] = str(source.with_suffix(".wl").absolute())
        execute.assert_called_once_with(str(self.kernel), [str(self.kernel), *expected])
        self.assertEqual(os.environ["RUSTFLOW_WOLFRAM_SCRIPT"], str(self.script))
        self.assertEqual(json.loads(os.environ["RUSTFLOW_WOLFRAM_SCRIPT_ARGUMENTS"]), arguments[4:])
        self.assertEqual(arguments[3], str(self.script), "the caller's list is unchanged")

    def test_child_invocation_replaces_only_its_own_script_environment(self):
        child = self.directory / "child.wl"
        child.write_text("Exit[0];\n")
        with patch.object(launcher.os, "execv"):
            launcher.launch(["-script", str(self.script), "parent"])
            launcher.launch(["-script", str(child), "child"])
        self.assertEqual(os.environ["RUSTFLOW_WOLFRAM_SCRIPT"], str(child))
        self.assertEqual(json.loads(os.environ["RUSTFLOW_WOLFRAM_SCRIPT_ARGUMENTS"]), ["child"])
        self.assertEqual(os.environ["RUSTFLOW_WOLFRAM_KERNEL"], str(self.kernel))

    def test_licensed_launcher_symlink_identity_is_preserved(self):
        link = self.directory / "chosen launcher"
        link.symlink_to(self.kernel)
        os.environ["RUSTFLOW_WOLFRAM_KERNEL"] = str(link)
        with patch.object(launcher.os, "execv") as execute:
            launcher.launch(["-script", str(self.script)])
        self.assertEqual(execute.call_args.args[0], str(link))
        self.assertEqual(execute.call_args.args[1][0], str(link))

    def test_missing_script_kernel_and_recursive_wrapper_are_rejected(self):
        for arguments in [[], ["-script"], ["-script", str(self.directory / "missing.wl")]]:
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                launcher.launch(arguments)
        del os.environ["RUSTFLOW_WOLFRAM_KERNEL"]
        with self.assertRaises(ValueError):
            launcher.launch(["-script", str(self.script)])
        # The entrypoint is normally executable after installation/copying.
        os.environ["RUSTFLOW_WOLFRAM_KERNEL"] = str(source)
        with patch.object(launcher.os, "access", return_value=True), self.assertRaisesRegex(ValueError, "this wrapper"):
            launcher.launch(["-script", str(self.script)])

    def test_nonexecutable_kernel_is_rejected_before_exec(self):
        self.kernel.chmod(0o644)
        with patch.object(launcher.os, "execv") as execute, self.assertRaises(ValueError):
            launcher.launch(["-script", str(self.script)])
        execute.assert_not_called()


if __name__ == "__main__":
    unittest.main()
