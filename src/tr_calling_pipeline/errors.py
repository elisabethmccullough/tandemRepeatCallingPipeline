"""Typed pipeline failures with machine-readable details."""

class PipelineError(RuntimeError):
    def __init__(self, message: str, **details: object):
        super().__init__(message)
        self.details = details

class ConfigurationError(PipelineError): pass
class ToolResolutionError(PipelineError): pass
class ToolVersionError(PipelineError): pass
class CommandExecutionError(PipelineError): pass
class CommandTimeoutError(CommandExecutionError): pass
class ExpectedOutputMissingError(CommandExecutionError): pass
class OutputConflictError(PipelineError): pass
class StageResumeError(PipelineError): pass
class ProvenanceError(PipelineError): pass
