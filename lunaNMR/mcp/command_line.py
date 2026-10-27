# ABOUTME: The `lunanmr-mcp` command-line parser, kept free of the MCP SDK so tests and docs can read it.
# ABOUTME: Its flags are what docs/MCP.md may mention.

import argparse


def build_parser():
    parser = argparse.ArgumentParser(
        prog='lunanmr-mcp', description='MCP server for lunaNMR, confined to one data folder')
    parser.add_argument('--root',
                        help='Dataset folder; the only place the server reads or writes')
    parser.add_argument('--check', action='store_true',
                        help='Verify this installation can run the server, then exit')
    return parser
