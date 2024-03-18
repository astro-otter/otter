'''
Custom exceptions for otter
'''

class FailedQuery(ValueError):
    def __str__(self):
        return "You're query/search did not return any results! Try again with different parameters!"

class IOError(ValueError):
    pass

class OtterLimitation(Exception):
    def __init__(self, msg):
        self.msg = "Current Limitation Found: " + msg

    def __str__(self):
        return self.msg

def TransientMergeError(Exception):
    pass
