# ABOUTME: Shared helper for the MCP server tests: one tool call through the SDK's in-process client.
# ABOUTME: Returns the tool's JSON envelope exactly as a client would receive it.

import asyncio
import json

from mcp.client import Client

from lunaNMR.mcp.data_root import DataRoot
from lunaNMR.mcp.server import build_server


def call_result(root, tool, arguments):
    """The protocol-level result: content, structured content and the error flag."""
    async def run():
        async with Client(build_server(DataRoot(root))) as client:
            return await client.call_tool(tool, arguments)
    return asyncio.run(run())


def call_tool(root, tool, arguments):
    return json.loads(call_result(root, tool, arguments).content[0].text)


def list_tools(root):
    async def run():
        async with Client(build_server(DataRoot(root))) as client:
            return (await client.list_tools()).tools
    return asyncio.run(run())
