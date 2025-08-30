# custom NoModulesException
class NoModulesException(Exception):
    def __init__(self, message):
        super().__init__(message)