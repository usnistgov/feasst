*************************************
Server
*************************************

Run FEASST as a server.
Note that the ports may not be closed and opened nicely.
Sometimes, the ports get a little fussy for some time and returns "ConnectionRefusedError: [Errno 111] Connection refused".
This error often seems to occur when a port does not close properly and waiting a few seconds for it to clean up will make the issue go away.
Or use a different port.

.. toctree::
   :glob:

   tutorial/tutorial*

FEASST plugin dependencies
============================

* monte_carlo

API
===

.. toctree::

   doc/toc
