import unittest
from util import BurnManTest
import tempfile
import os

from burnman.utils.misc import extract_lines_between_markers
from burnman.utils.misc import run_cli_program_with_input


class test_utils(BurnManTest):
    def test_run_cli(self):
        output = run_cli_program_with_input(["cat"], "Hello, CLI!\n", verbose=False)
        self.assertEqual(output, "Hello, CLI!\n")

    def test_extract_lines(self):
        test_content = """Header
    START_MARK
    Line 1
    Line 2
    END_MARK
    Footer"""

        # Create a temporary file with the test content
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp:
            tmp.write(test_content)
            tmp_path = tmp.name

        # Exclusive test (default)
        result = extract_lines_between_markers(tmp_path, "START_MARK", "END_MARK")
        self.assertEqual(result[0], "Line 1")
        self.assertEqual(result[1], "Line 2")

        # Inclusive test
        result = extract_lines_between_markers(
            tmp_path, "START_MARK", "END_MARK", inclusive=True
        )
        self.assertEqual(result[0], "START_MARK")
        self.assertEqual(result[1], "Line 1")
        self.assertEqual(result[2], "Line 2")
        self.assertEqual(result[3], "END_MARK")

        os.remove(tmp_path)  # Clean up


if __name__ == "__main__":
    unittest.main()
