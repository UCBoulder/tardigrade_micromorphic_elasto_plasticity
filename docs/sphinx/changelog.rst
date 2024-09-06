.. _changelog:


#########
Changelog
#########

******************
1.3.3 (unreleased)
******************

Bug Fixes
=========
- Enabled use of new micromorphic_tools library (:pull:`8`). By `Nathan Miller`_.
- Corrected missing catch for convergence errors for the non-Jacobian hydra-model version (:pull:`10`). By `Nathan Miller`_.

Internal Changes
================
- Updated to use new error_tools nested exception message catch utility (:pull:`10`). By `Nathan Miller`_.
- Added flag to allow for all errors to be treated as convergence errors (:pull:`11`). By `Nathan Miller`_.
- Enable the Hydra pre-conditioner (:pull:`12`). By `Nathan Miller`_.
- Accommodate new Hydra base class inputs (:pull:`14`). By `Nathan Miller`_.
- Accommodate new integration scheme for the plasticity (:pull:`15`). By `Nathan Miller`_.
- Rolled back new integration scheme (:pull:`16`). By `Nathan Miller`_.
- Removed trailing whitespace from add_library in CMakeLists file (:pull:`17`). By `Nathan Miller`_.
- Added ability to set package versions when doing FetchContent builds (:pull:`18`). By `Nathan Miller`_.
- Enabled correct functionality with new Levenberg-Marquardt backup solver (:pull:`19`). By `Nathan Miller`_.
- Added tests of difficult problems (:pull:`19`). By `Nathan Miller`_.
- Added optimization-based solver for difficult cases (:pull:`20`). By `Nathan Miller`_.

New Features
============
- Enabled timestep cutback requests from hydra (:pull:`5`). By `Nathan Miller`_.
- Set the plastic deformation to be fully implicitly integrated (:pull:`6`). By `Nathan Miller`_.
- Added additional information to output error messages (:pull:`13`). By `Nathan Miller`_.

******************
1.3.2 (2024-01-24)
******************

Release
=======
- Released version (:pull:`4`). By `Nathan Miller`_.

New Features
============
- Added a tardigrade_hydra based version of the elasto-plastic model (:pull:`3`). By `Nathan Miller`_.

Internal Changes
================
- Enabled use of github actions to run tests (:pull:`1`). By `Nathan Miller`_.
- Added tardigrade hydra as a dependency (:pull:`2`). By `Nathan Miller`_.
- Updated tests to work with hydra 0.4.x (:pull:`4`). By `Nathan Miller`_.

******************
1.3.1 (2023-07-25)
******************

Breaking Changes
================
- Change project, package, and namespace to use the 'tardigrade' prefix (:issue:`6`, :merge:`13`). By `Kyle Brindley`_.

******************
1.2.1 (2023-07-12)
******************

Internal Changes
================
- Replace build scripts with direct use of CMake commands in CI configuration (:issue:`2`, :merge:`8`). By `Kyle
  Brindley`_.
- Create CI environment (:issue:`3`, :merge:`9`). By `Kyle Brindley`_.
- Use setuptools_scm for version number (:issue:`4`, :merge:`10`). By `Kyle Brindley`_.
- Add conda package recipe and deploy CI jobs (:issue:`5`, :merge:`11`). By `Kyle Brindley`_.

******************
1.1.0 (2022-08-16)
******************

- Moved the code to the cpp_stub format (:merge:`1`). By `Nathan Miller`_.
- Moved the tests to the BOOST test format (:merge:`2`). By `Nathan Miller`_.
- Removed old material library interface definitions (:merge:`3`). By `Nathan Miller`_.
