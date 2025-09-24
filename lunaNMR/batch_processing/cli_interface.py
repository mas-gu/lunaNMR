#!/usr/bin/env python3
"""
Command Line Interface for LunaNMR Batch Processing

This module provides a comprehensive CLI interface for batch processing NMR spectra.
It supports configuration files, preset selection, parameter optimization, and
detailed progress reporting.

Features:
- Complete argument parsing and validation
- Support for configuration files and presets
- Real-time progress reporting
- Flexible output options
- Error handling and recovery options
"""

import argparse
import sys
import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
import json
import signal

# Import batch processing components
try:
    from .batch_processor import BatchProcessor
    from .parameter_optimizer import ParameterOptimizer
    from .config_manager import ConfigManager
    BATCH_COMPONENTS_AVAILABLE = True
except ImportError:
    BATCH_COMPONENTS_AVAILABLE = False

class CLIInterface:
    """
    Command-line interface for batch NMR processing.

    Provides a user-friendly interface for processing multiple NMR spectra
    with automatic parameter optimization and flexible configuration options.
    """

    def __init__(self):
        """Initialize the CLI interface."""
        self.logger = logging.getLogger(__name__)
        self.interrupted = False

        # Set up signal handling for graceful shutdown
        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)

    def _signal_handler(self, signum, frame):
        """Handle interrupt signals for graceful shutdown."""
        print("\n\n⚠️  Processing interrupted by user")
        print("Gracefully shutting down...")
        self.interrupted = True
        sys.exit(0)

    def create_parser(self) -> argparse.ArgumentParser:
        """
        Create the argument parser for the CLI interface.

        Returns:
            Configured ArgumentParser instance
        """
        parser = argparse.ArgumentParser(
            description="LunaNMR Batch Processing - Automated NMR Spectrum Analysis",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Examples:
  # Basic usage - process all spectra in folder
  python -m lunaNMR.batch_processing /path/to/spectra

  # Specify nucleus type for all spectra
  python -m lunaNMR.batch_processing /path/to/spectra --nucleus 15N1H

  # Use a preset configuration
  python -m lunaNMR.batch_processing /path/to/spectra --preset conservative

  # Use custom configuration file
  python -m lunaNMR.batch_processing /path/to/spectra --config config.json

  # Enable parameter optimization
  python -m lunaNMR.batch_processing /path/to/spectra --optimize

  # Generate example configuration files
  python -m lunaNMR.batch_processing --create-examples /path/to/output
            """
        )

        # Primary arguments
        parser.add_argument(
            'folder',
            nargs='?',
            help='Path to folder containing NMR spectra (required unless using --create-examples)'
        )

        # Processing options
        parser.add_argument(
            '--nucleus',
            choices=['15N1H', '13C1H', '1H'],
            help='Fixed nucleus type (if not specified, will auto-detect from filenames)'
        )

        parser.add_argument(
            '--preset',
            choices=['15N1H', '13C1H', '1H', 'conservative', 'aggressive', 'high_throughput'],
            help='Use a preset configuration optimized for specific conditions'
        )

        parser.add_argument(
            '--config',
            type=str,
            metavar='FILE',
            help='Path to custom JSON configuration file'
        )

        parser.add_argument(
            '--optimize',
            action='store_true',
            help='Enable automatic parameter optimization for each spectrum'
        )

        # Output and logging options
        parser.add_argument(
            '--output',
            type=str,
            metavar='DIR',
            help='Directory for output files (default: current directory)'
        )

        parser.add_argument(
            '--log-level',
            choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
            default='INFO',
            help='Set logging level (default: INFO)'
        )

        parser.add_argument(
            '--quiet',
            action='store_true',
            help='Minimal console output (only errors and final summary)'
        )

        parser.add_argument(
            '--verbose',
            action='store_true',
            help='Verbose output with detailed progress information'
        )

        # File handling options
        parser.add_argument(
            '--extensions',
            nargs='+',
            default=['.ft2', '.ft', '.pipe', '.ucsf', '.nmrpipe'],
            help='File extensions to process (default: .ft2 .ft .pipe .ucsf .nmrpipe)'
        )

        parser.add_argument(
            '--recursive',
            action='store_true',
            help='Search for files recursively in subdirectories'
        )

        # Error handling options
        parser.add_argument(
            '--skip-errors',
            action='store_true',
            help='Continue processing even if individual spectra fail (default behavior)'
        )

        parser.add_argument(
            '--strict',
            action='store_true',
            help='Stop processing on first error'
        )

        # Utility options
        parser.add_argument(
            '--create-examples',
            type=str,
            metavar='DIR',
            help='Create example configuration files in specified directory'
        )

        parser.add_argument(
            '--list-presets',
            action='store_true',
            help='List available preset configurations and exit'
        )

        parser.add_argument(
            '--dry-run',
            action='store_true',
            help='Show what would be processed without actually processing'
        )

        parser.add_argument(
            '--version',
            action='version',
            version='LunaNMR Batch Processing v1.0.0'
        )

        return parser

    def validate_arguments(self, args: argparse.Namespace) -> bool:
        """
        Validate command line arguments.

        Args:
            args: Parsed command line arguments

        Returns:
            True if arguments are valid, False otherwise
        """
        # Check for conflicting options
        if args.quiet and args.verbose:
            print("Error: --quiet and --verbose options are mutually exclusive")
            return False

        if args.skip_errors and args.strict:
            print("Error: --skip-errors and --strict options are mutually exclusive")
            return False

        if args.preset and args.config:
            print("Error: --preset and --config options are mutually exclusive")
            return False

        # Check for required folder argument
        if not args.create_examples and not args.list_presets and not args.folder:
            print("Error: folder argument is required unless using --create-examples or --list-presets")
            return False

        # Check folder exists
        if args.folder and not Path(args.folder).exists():
            print(f"Error: Folder '{args.folder}' does not exist")
            return False

        # Check config file exists
        if args.config and not Path(args.config).exists():
            print(f"Error: Configuration file '{args.config}' does not exist")
            return False

        return True

    def setup_logging(self, args: argparse.Namespace):
        """
        Set up logging based on command line arguments.

        Args:
            args: Parsed command line arguments
        """
        # Determine log level
        if args.quiet:
            console_level = logging.ERROR
        elif args.verbose:
            console_level = logging.DEBUG
        else:
            console_level = getattr(logging, args.log_level)

        # Configure root logger
        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[]
        )

        # Console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(console_level)

        if args.quiet:
            console_format = '%(message)s'
        elif args.verbose:
            console_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        else:
            console_format = '%(levelname)s - %(message)s'

        console_handler.setFormatter(logging.Formatter(console_format))

        # Add handlers
        logger = logging.getLogger()
        logger.addHandler(console_handler)

    def load_configuration(self, args: argparse.Namespace) -> Dict[str, Any]:
        """
        Load configuration from arguments, preset, or config file.

        Args:
            args: Parsed command line arguments

        Returns:
            Configuration dictionary
        """
        config_manager = ConfigManager()

        # Start with default configuration
        if args.preset:
            config = config_manager.get_preset_config(args.preset)
            print(f"Using preset configuration: {args.preset}")
        elif args.config:
            config = config_manager.load_config(args.config)
            print(f"Using configuration file: {args.config}")
        else:
            config = config_manager.default_config.copy()
            print("Using default configuration")

        # Override with command line arguments
        if args.extensions:
            config['file_handling']['extensions'] = args.extensions

        if args.recursive:
            config['file_handling']['recursive_search'] = True

        if args.skip_errors:
            config['processing']['skip_on_error'] = True
        elif args.strict:
            config['processing']['skip_on_error'] = False

        if args.optimize:
            config['optimization']['enable_auto_optimization'] = True

        # Logging configuration
        config['logging']['level'] = args.log_level
        config['logging']['detailed_logging'] = args.verbose

        return config

    def create_example_files(self, output_dir: str):
        """
        Create example configuration files.

        Args:
            output_dir: Directory to create example files in
        """
        config_manager = ConfigManager()
        output_path = Path(output_dir)

        print(f"Creating example configuration files in: {output_path}")

        # Create preset configuration files
        created_files = config_manager.create_preset_config_files(output_path)

        # Create main example configuration
        example_file = output_path / "example_config.json"
        config_manager.create_example_config(example_file)
        created_files.append(example_file)

        print(f"Created {len(created_files)} configuration files:")
        for file_path in created_files:
            print(f"  - {file_path.name}")

        print("\nYou can now use these configurations with:")
        print("  python -m lunaNMR.batch_processing /path/to/spectra --config example_config.json")
        print("  python -m lunaNMR.batch_processing /path/to/spectra --preset 15N1H")

    def list_presets(self):
        """List available preset configurations."""
        config_manager = ConfigManager()
        presets = config_manager.get_available_presets()

        print("Available preset configurations:\n")
        for name, description in presets.items():
            print(f"  {name:15} - {description}")

        print("\nUsage:")
        print("  python -m lunaNMR.batch_processing /path/to/spectra --preset PRESET_NAME")

    def perform_dry_run(self, args: argparse.Namespace, config: Dict[str, Any]):
        """
        Perform a dry run showing what would be processed.

        Args:
            args: Parsed command line arguments
            config: Configuration dictionary
        """
        print("=== DRY RUN - No actual processing will be performed ===\n")

        # Create temporary processor to find files
        processor = BatchProcessor(config)

        try:
            spectrum_files = processor.find_spectrum_files(args.folder)

            print(f"Found {len(spectrum_files)} spectrum files:")
            for i, file_path in enumerate(spectrum_files[:10]):  # Show first 10
                nucleus_type = args.nucleus or processor.detect_nucleus_type(file_path.name)
                print(f"  {i+1:3d}. {file_path.name} (detected: {nucleus_type})")

            if len(spectrum_files) > 10:
                print(f"  ... and {len(spectrum_files) - 10} more files")

            print(f"\nConfiguration:")
            print(f"  Nucleus type: {args.nucleus or 'auto-detect'}")
            print(f"  Optimization: {'enabled' if args.optimize else 'disabled'}")
            print(f"  Error handling: {'skip errors' if config['processing']['skip_on_error'] else 'strict'}")
            print(f"  File extensions: {', '.join(config['file_handling']['extensions'])}")

        except Exception as e:
            print(f"Error during dry run: {e}")

    def run_batch_processing(self, args: argparse.Namespace, config: Dict[str, Any]) -> bool:
        """
        Run the actual batch processing.

        Args:
            args: Parsed command line arguments
            config: Configuration dictionary

        Returns:
            True if processing was successful, False otherwise
        """
        try:
            # Create processor
            processor = BatchProcessor(config)

            # Add parameter optimizer if requested
            if args.optimize:
                optimizer = ParameterOptimizer(config.get('optimization', {}))
                # The processor would use the optimizer internally
                print("Parameter optimization enabled")

            # Progress reporting
            if not args.quiet:
                print(f"Starting batch processing of: {args.folder}")
                print(f"Target nucleus type: {args.nucleus or 'auto-detect'}")

            # Process the folder
            results = processor.process_folder(
                folder_path=args.folder,
                nucleus_type=args.nucleus,
                auto_optimize=args.optimize
            )

            # Report final results
            if not args.quiet:
                print("\n" + "="*60)
                print("BATCH PROCESSING COMPLETED")
                print("="*60)
                print(f"Total files processed: {results['processed_files']}/{results['total_files']}")
                print(f"Success rate: {(results['processed_files']/results['total_files']*100):.1f}%")
                print(f"ML training samples collected: {results['total_ml_samples']}")

                if results['failed_files'] > 0:
                    print(f"Failed files: {results['failed_files']}")

            return results['processed_files'] > 0

        except KeyboardInterrupt:
            print("\nProcessing interrupted by user")
            return False
        except Exception as e:
            print(f"Error during batch processing: {e}")
            return False

    def main(self, argv: Optional[List[str]] = None) -> int:
        """
        Main CLI entry point.

        Args:
            argv: Optional command line arguments (uses sys.argv if None)

        Returns:
            Exit code (0 for success, 1 for error)
        """
        if not BATCH_COMPONENTS_AVAILABLE:
            print("Error: Batch processing components not available")
            print("Please ensure lunaNMR core components are properly installed")
            return 1

        # Parse arguments
        parser = self.create_parser()
        args = parser.parse_args(argv)

        # Validate arguments
        if not self.validate_arguments(args):
            return 1

        # Set up logging
        self.setup_logging(args)

        try:
            # Handle utility commands
            if args.create_examples:
                self.create_example_files(args.create_examples)
                return 0

            if args.list_presets:
                self.list_presets()
                return 0

            # Load configuration
            config = self.load_configuration(args)

            # Handle dry run
            if args.dry_run:
                self.perform_dry_run(args, config)
                return 0

            # Run batch processing
            success = self.run_batch_processing(args, config)
            return 0 if success else 1

        except Exception as e:
            print(f"Unexpected error: {e}")
            return 1

def main():
    """Entry point for command line usage."""
    cli = CLIInterface()
    sys.exit(cli.main())

if __name__ == "__main__":
    main()