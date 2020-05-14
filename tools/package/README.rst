Package scripts
**************************

The package scripts help a user to allow certain files to diverge from or stay hidden from the central FEASST repository.
For example, perhaps a newly implemented method is unpublished and not ready for public distribution.
One option is to simply make a private branch within your own private repository that you will not push to the public repository.
This method is fine if you only want to receive updates.
But if you want to share some of the development within your private branch, but not all, this becomes more complicated.

The other option is to create separate "public" and "private" repositories.
The example package directory `<tools/package>`_ contains scripts and example files to simplify the synchronization of these two separate repositories.

* `package.sh`_: send updates in private repository to the public repository, after filtering.
* `unpack.sh`_: receive changes from public repository into private repository, after filtering.
* `private_files.txt`_: list of files or directories that will not be shared to, or altered by, the public repository.
* `stubs`_: directory of files, which will overwrite current files when sent to public repository, and will not be changed by unpack.
* `feasst_private`_: a directory with a file of this name will behave as if listed in `private_files.txt`_.

package.sh
=============

.. code-block:: bash

   ./package.sh /path/to/public/repo

This script sends the contents of the current FEASST repository to a different FEASST repository, after filtering of files based on `private_files.txt` and `stubs`.
At the end, it lists the differences between the two repositories.
If, for example, a file is deleted or renamed in your private repository and sent to the public one, this will show up as a difference in the repositories.

unpack.sh
===========

.. code-block:: bash

   ./unpack.sh /path/to/public/repo

This script receives the contents of a different FEASST repository, after filtering of files based on `private_files.txt` and `stubs`.

private_files.txt
==================

List each file or directory that you do not want `package.sh` to send to another repository.
In addition, `unpack.sh` will not receive these files from the public repository.
Finally, you may also add an empty file, `feasst_private`_ to any directory, and the pack and unpack scripts will append these directories onto the private_files list.

stubs
======

In the stubs directory, include files which you would like to replace the existing file in your repository.
For example, you may select some custom options for CMakeLists.txt but not want to push those back to the public repository.
Thus, the default CMakeLists.txt may be included as a stub which is copied on top of your own implementation when sending changes.

feasst_private
===============

Create a file with this name in order to keep the contents of the entire directory private.
`package.sh` will behave as if this directory was listed in `private_files.txt`_.


