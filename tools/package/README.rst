Package scripts
**************************

The package scripts help a user to allow certain files to diverge from or stay hidden from the central FEASST respository.
For example, perhaps a newly implemented method is unpublished and not ready for public distribution.
One option is to simply make a private branch within your own private repository that you will not push to the public respository.
This method is fine if you only which to receive updates.
But if you want to share some of the development within your private branch, but not all, this becomes more complicated.

The other option is to create separate "public" and "private" respositories.
The example package directory `<tools/package>`_ contains scripts and example files to simplify the synchronization of these two separate respositories.

* `package.sh`_: send updates in private respository to the public respository, after filtering.
* `unpack.sh`_: receive changes from public respository into private respository, after filtering.
* `private_files.txt`_: list of files or directories that will not be shared to, or altered by, the public respository.
* `stubs`_: directory of files, which will overwrite current files when sent to public respository, and will not be changed by unpack.

package.sh
=============

.. code-block:: bash

   ./package.sh /path/to/public/repo

This script sends the contents of the current FEASST repository to a different FEASST repository, after filtering of files based on `private_files.txt` and `stubs`.
At the end, it lists the differences between the two respositories.
If, for example, a file is deleted or renamed in your private respository and sent to the public one, this will show up as a difference in the respositories.

unpack.sh
===========

.. code-block:: bash

   ./unpack.sh /path/to/public/repo

This script receives the contents of a different FEASST respository, after filtering of files based on `private_files.txt` and `stubs`.

private_files.txt
==================

List each file or directory that you do not want `package.sh` to send to another respository.
In addition, `unpack.sh` will not receive these files from the public repository.

stubs
======

In the stubs directory, include files which you would like to replace the existing file in your repository.
For example, you may select some custom options for CMakeLists.txt but not want to push those back to the public repository.
Thus, the default CMakeLists.txt may be included as a stub which is copied on top of your own implementation when sending changes.


