***********************************************
Modify FEASST
***********************************************

This plugin contains example templates for creating your own custom FEASST Model, Action or Analyze.
Users can make their own class by copy-pasting and renaming a pair of header (.h) and implementation (.cpp) files.
Find an existing class that does something as close as possible to what you want to do, copy and rename it (as detailed below), make some changes and test it.
If it turns out to be useful, you can become a developer/contributor by submitting a pull request on GitHub.

Specifically, these are the steps required to copy and rename an existing class in FEASST.

1. Copy and paste this h file, and its corresponding cpp file with a new name. Replace ModelExample with NewName in the cpp and h files, and MODEL_EXAMPLE with NEW_NAME in the h file.

.. code-block:: bash

    cd /path/to/feasst/plugin/example
    sed "s/MODEL_EXAMPLE/NEW_NAME/g" include/model_example.h > include/new_name.h
    sed "s/model_example\.h/new_name.h/g" src/model_example.cpp > src/new_name.cpp
    sed -i "s/ModelExample/NewName/g" include/new_name.h src/new_name.cpp

2. Reinstall FEASST, starting with the cmake command.

3. (optional) Add C++ tests by performing steps 1 and 2 for /path/to/feasst/plugin/example/test/model_example.cpp . Add tutorials in /path/to/feasst/plugin/example/tutorial .

4. (optional) Change the serialization version found in two places in the cpp
   file from the current 5023 to something else (in both places).

This process can be used on any class in FEASST, not just the ones in the example plugin.
Get started by reading /path/to/feasst/plugin/example/model_example.h

In addition, use the example plugin as a template to create new plugins with the following steps:

1. Copy and rename the entire plugin.

.. code-block:: bash

    cd /path/to/feasst/plugin
    cp -r example new_plugin

2. Remove the classes (.h and .cpp files) that you do not need and make sure the ones you need are renamed to something unique, as described in the previous steps 1-3 for making a new class detailed above.

3. Find "set(FEASST_PLUGINS ...)" in /path/to/feasst/CMakeLists.txt and add "new_plugin" to the semicolon-separated list.

4. Add one class present in the new_plugin near the bottom of /path/to/feasst/py/depend.py. For example, if there is a class NewName in new_plugin, add the following line near all of the similar lines near the end of depend.py: "if 'new_plugin' in include_plugin: select_classes.append("NewName")

5. Reinstall FEASST, starting with the cmake command.

.. toctree::
   :glob:

   tutorial/tutorial*

FEASST plugin dependencies
============================

* system

API
===

.. toctree::

   doc/toc
