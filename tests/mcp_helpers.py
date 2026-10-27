# ABOUTME: Shared helper for the MCP server tests: one tool call through the SDK's in-process client.
# ABOUTME: Returns the tool's JSON envelope exactly as a client would receive it.

import asyncio
import json

from mcp.client import Client

from lunaNMR.mcp.data_root import DataRoot
from lunaNMR.mcp.server import build_server


def call_tool(root, tool, arguments):
    async def run():
        async with Client(build_server(DataRoot(root))) as client:
            result = await client.call_tool(tool, arguments)
            return json.loads(result.content[0].text)
    return asyncio.run(run())
