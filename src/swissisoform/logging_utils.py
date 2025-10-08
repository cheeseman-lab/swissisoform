"""Logging utilities for SwissIsoform.

This module provides utilities for structured logging with proper formatting
and context management, replacing hardcoded indentation and box-drawing characters.
"""

import logging
from contextlib import contextmanager
from typing import Optional


class DebugLogger:
    """Wrapper for logger with context-aware indentation.

    Provides clean, professional debug output without hardcoded
    indentation or box-drawing characters.

    Example:
        >>> debug_log = DebugLogger(logger, debug=True)
        >>> debug_log.info("Processing gene X")
        >>> with debug_log.section("Validating mutation"):
        ...     debug_log.debug("Genomic position: 12345")
        ...     debug_log.debug("Reference allele: A")
        ...     debug_log.success("Validation passed")
    """

    def __init__(self, logger: logging.Logger, debug: bool = False):
        """Initialize the debug logger.

        Args:
            logger: The underlying logger instance
            debug: Whether debug mode is enabled
        """
        self.logger = logger
        self._debug_enabled = debug
        self._indent_level = 0
        self._indent_str = "  "  # 2 spaces per level

    def _format_message(self, message: str, prefix: str = "") -> str:
        """Format a message with current indentation.

        Args:
            message: The message to format
            prefix: Optional prefix (e.g., "✓", "✗", "→")

        Returns:
            Formatted message string
        """
        indent = self._indent_str * self._indent_level
        if prefix:
            return f"{indent}{prefix} {message}"
        return f"{indent}{message}"

    @contextmanager
    def section(self, title: str):
        """Create an indented section for related log messages.

        Args:
            title: Section title

        Yields:
            None
        """
        if self._debug_enabled:
            self.logger.debug(self._format_message(title))
        self._indent_level += 1
        try:
            yield
        finally:
            self._indent_level -= 1

    def debug(self, message: str, prefix: str = ""):
        """Log a debug message if debug mode is enabled.

        Args:
            message: Message to log
            prefix: Optional prefix character/emoji
        """
        if self._debug_enabled:
            self.logger.debug(self._format_message(message, prefix))

    def info(self, message: str, prefix: str = ""):
        """Log an info message.

        Args:
            message: Message to log
            prefix: Optional prefix character/emoji
        """
        self.logger.info(self._format_message(message, prefix))

    def warning(self, message: str, prefix: str = "⚠"):
        """Log a warning message.

        Args:
            message: Message to log
            prefix: Optional prefix character/emoji (default: warning symbol)
        """
        self.logger.warning(self._format_message(message, prefix))

    def error(self, message: str, prefix: str = "✗"):
        """Log an error message.

        Args:
            message: Message to log
            prefix: Optional prefix character/emoji (default: X mark)
        """
        self.logger.error(self._format_message(message, prefix))

    def success(self, message: str, prefix: str = "✓"):
        """Log a success message at debug level.

        Args:
            message: Message to log
            prefix: Optional prefix character/emoji (default: check mark)
        """
        if self._debug_enabled:
            self.logger.debug(self._format_message(message, prefix))

    def field(self, name: str, value: any):
        """Log a field/value pair at debug level.

        Args:
            name: Field name
            value: Field value
        """
        if self._debug_enabled:
            self.logger.debug(self._format_message(f"{name}: {value}"))


def configure_logging(
    level: int = logging.INFO,
    format_string: Optional[str] = None,
    include_timestamp: bool = True,
) -> None:
    """Configure logging for SwissIsoform applications.

    Args:
        level: Logging level (e.g., logging.DEBUG, logging.INFO)
        format_string: Custom format string (if None, uses default)
        include_timestamp: Whether to include timestamps in logs
    """
    if format_string is None:
        if include_timestamp:
            format_string = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        else:
            format_string = "%(name)s - %(levelname)s - %(message)s"

    logging.basicConfig(
        level=level,
        format=format_string,
        datefmt="%Y-%m-%d %H:%M:%S",
    )


# Example usage in scripts:
if __name__ == "__main__":
    # Configure logging
    configure_logging(level=logging.DEBUG, include_timestamp=False)

    # Create a debug logger
    logger = logging.getLogger(__name__)
    debug_log = DebugLogger(logger, debug=True)

    # Example output
    debug_log.info("Processing gene TP53")

    with debug_log.section("Validating mutation"):
        debug_log.field("Genomic position", "chr17:7674220")
        debug_log.field("Reference allele", "C")
        debug_log.field("Alternate allele", "T")
        debug_log.success("Reference allele validated")

    with debug_log.section("Applying mutation"):
        debug_log.debug("Extracting coding sequence")
        debug_log.field("CDS length", "1182 bp")
        debug_log.field("Protein length", "393 AA")
        debug_log.success("Mutation applied successfully")

    debug_log.info("Processing complete")
