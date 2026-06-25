import threading

class ExcThread(threading.Thread):
    """
    Handle exceptions in threads.

    >>> from feasst import excthread
    >>> def failing_task(): raise RuntimeError("Error")
    >>> thread = ExcThread(target=failing_task)
    >>> thread.start()
    >>> thread.join()
    >>> if thread.exception: print(thread.exception)
    Error
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.exception = None

    def run(self):
        try:
            if self._target:
                self._target(*self._args, **self._kwargs)
        except Exception as e:
            self.exception = e

if __name__ == "__main__":
    import doctest
    doctest.testmod()
