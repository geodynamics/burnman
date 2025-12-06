2021-11-02
^^^^^^^^^^

* Running the examples now requires BurnMan to be installed. This change is
  made to ensure that tester behaviour always matches installed behaviour,
  rather than local behaviour (the cause of the problem with the 1.0.0 release).
  If the user wishes to develop the BurnMan code, they may install the module
  from the top-level directory using pip with the `-e` (editable) flag:
  `python -m pip install -e .`.

  *Bob Myhill, 2021/11/02*
