# Contributing to BurnMan
BurnMan is a community project that lives by the participation of its
members — i.e., including you! It is our goal to build an inclusive
and participatory community so we are happy that you are interested in
participating! 

## Getting started with git and GitHub
GitHub provides a helpful
guide on the process of contributing to an open-source project
[here](https://opensource.guide/how-to-contribute/).

## Asking and answering questions about BurnMan
For questions about BurnMan on all levels, please use the 
[BurnMan forum](https://community.geodynamics.org/c/burnman).

## Bug reports
It is a great help to the community if you report any bugs that you
may find. We keep track of all open issues related to BurnMan
[here](https://github.com/geodynamics/burnman/issues). 

Please follow these simple instructions before opening a new bug report:

- Do a quick search in the list of open and closed issues for a duplicate of
  your issue.
- If you did not find an answer, open a new
  [issue](https://github.com/geodynamics/burnman/issues/new) and explain your
  problem in as much detail as possible.
- Attach as much as possible of the following information to your issue:
  - a minimal script that reproduces the issue,
  - the error message you saw on your screen,
  - any information that helps us understand why you think this is a bug, and
    how to reproduce it.

## Making BurnMan better
BurnMan is a community project, and we are encouraging all kinds of
contributions. Obvious candidates are bugfixes and implementations of new datasets or
fitting routines. Other much appreciated contributions are new examples,
tests, benchmarks, fixing typos or updating outdated documentation. 
If you consider making a large addition or change to core functionality,
please open a new [issue](https://github.com/geodynamics/burnman/issues/new)
first, to discuss your idea with one of the maintainers. This allows us to give you early
feedback and prevents you from spending much time on a project that might already be
planned, or that conflicts with other plans.

### Opening pull requests

To make a change to BurnMan you should:
- Create a
[fork](https://guides.github.com/activities/forking/#fork) (through GitHub) of
the code base.
- Create a new
[branch](https://guides.github.com/introduction/flow/) (sometimes called a
feature branch) on which to make your modifications.
- You can propose that your branch be merged into the BurnMan
code by opening a [pull request](https://guides.github.com/introduction/flow/).
This will give others a chance to review your code. 

We follow the philosophy that no pull request (independent of the author) is
merged without a review from one other member of the community, and approval of
one of the maintainers. This applies to maintainers as well as to first-time
contributors. We know that review can be a daunting process, but pledge to
keep all comments friendly and supportive. We are as
interested in making BurnMan better as you are!

While this seems very
formal, keeping all of the code review in one place makes it easier to
coordinate changes to the code. Please do
not hesitate to ask questions about the workflow on the forum if you are
not sure what to do.

### Coding conventions

Since BurnMan is a growing project with several contributors we
use [black](https://black.readthedocs.io/en/stable/)
to keep the style of the source code
consistent. See the [project file](pyproject.toml) for the most recent version.
Indentation and other code standardisation is achieved by running the
[indent script](contrib/utilities/indent). If you
are new to the project then we will work with you to ensure your contributions
are formatted with this style, so please do not think of it as a road block if
you would like to contribute some code.

### Changelog entries

If your new pull request creates a change that is noticeable to BurnMan users,
please add a new file to our [changelog](docs/changelog) folder.
The name of the file consists of the date of the change
(approximately) and the name of the author. Use the existing files for a
guide to the format. This ensures you will get credit for your work.

## Acknowledgment of contributions

The BurnMan community is grateful for every contribution! But, beyond
this gratitude, there are also several *formal*
ways in which your contribution will be acknowledged by the BurnMan community:
- Every commit that is merged into the BurnMan repository makes you part of
  the group of 
  [contributors](https://github.com/geodynamics/burnman/graphs/contributors).
- For every release the most significant entries of our
  [changelog](docs/changelog) are selected to generate our release announcements.
  Additionally, all entries remain available within the
  [docs/changelog](docs/changelog) directory.
- If you contributed a significant part of the manual, you will be listed as
  one of the contributing authors of the manual.
- Significant contributions will lead to your name being included in the
  [AUTHORS](AUTHORS.md) file in the main repository. Criteria for inclusion:

  - A profound understanding of BurnMan's structure and vision;
  - A proven willingness to further the project's goals and help other users;
  - Significant contributions to BurnMan (not necessarily only source code,
    also mailing list advice, documentation, benchmarks, tutorials);
  - Contributions to BurnMan for more than one year.

## License
BurnMan is published under the [GPL v2 or newer](LICENSE); while you
will retain copyright on your contributions, all changes to the code
must be provided under this common license.
