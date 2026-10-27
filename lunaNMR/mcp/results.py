# ABOUTME: The shape of every MCP tool result: the CLI envelope, an error flag, and a size bound.
# ABOUTME: Also the read/write annotations tools declare, so a client can tell what each one changes.

import json

from mcp.types import CallToolResult, TextContent, ToolAnnotations

MAX_RESULT_CHARS = 200_000


def result(envelope):
    """A tool result whose error flag follows the envelope's `ok`, so a client treats a
    failed run as a failure and still reads the reason. Anything larger than a client should
    be handed is replaced by an error: the text comes from files the dataset controls."""
    text = json.dumps(envelope)
    if len(text) > MAX_RESULT_CHARS:
        envelope = {'ok': False, 'error': {
            'type': 'ResultTooLarge',
            'message': f'the result was {len(text)} characters; read the files in the output '
                       f'folder instead'}}
        text = json.dumps(envelope)
    return CallToolResult(
        content=[TextContent(type='text', text=text)],
        structured_content=envelope,
        is_error=not envelope.get('ok', False))


def failure(error):
    return result({'ok': False, 'error': {'type': type(error).__name__, 'message': str(error)}})


READS_ONLY = ToolAnnotations(read_only_hint=True, open_world_hint=False)
# Rewrites the residue list and settings file beside its input (merging, so decisions survive).
REWRITES_BESIDE_INPUT = ToolAnnotations(
    read_only_hint=False, destructive_hint=True, idempotent_hint=True, open_world_hint=False)
CANCELS_A_JOB = ToolAnnotations(
    read_only_hint=False, destructive_hint=True, idempotent_hint=True, open_world_hint=False)
WRITES_NEW_FOLDER = ToolAnnotations(
    read_only_hint=False, destructive_hint=False, idempotent_hint=False, open_world_hint=False)
