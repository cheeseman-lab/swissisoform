"""Configuration management for SwissIsoform.

This module loads configuration from environment variables or .env file.
Sensitive credentials like API keys should be stored in .env (not committed to git).
"""

import os
from pathlib import Path
from typing import Optional


def load_env_file(env_path: Optional[str] = None) -> None:
    """Load environment variables from .env file.

    Args:
        env_path: Optional path to .env file. If not provided, searches in standard locations.
    """
    if env_path:
        env_file = Path(env_path)
    else:
        # Search for .env in project root (up to 3 levels from this file)
        current = Path(__file__).parent
        for _ in range(4):
            env_file = current / ".env"
            if env_file.exists():
                break
            current = current.parent
        else:
            # .env file not found, will use defaults
            return

    if not env_file.exists():
        return

    # Simple .env parser (doesn't require python-dotenv dependency)
    with open(env_file, "r") as f:
        for line in f:
            line = line.strip()
            # Skip comments and empty lines
            if not line or line.startswith("#"):
                continue

            # Parse KEY=VALUE
            if "=" in line:
                key, value = line.split("=", 1)
                key = key.strip()
                value = value.strip()
                # Don't override existing environment variables
                if key not in os.environ:
                    os.environ[key] = value


# Load .env file on module import
load_env_file()


class Config:
    """Configuration class for SwissIsoform API credentials and endpoints."""

    # NCBI E-utilities API Key
    # Get your free API key at: https://www.ncbi.nlm.nih.gov/account/
    NCBI_API_KEY: Optional[str] = os.getenv("NCBI_API_KEY")

    # gnomAD API endpoint
    GNOMAD_API_URL: str = os.getenv(
        "GNOMAD_API_URL", "https://gnomad.broadinstitute.org/api"
    )

    # ClinVar/NCBI E-utilities base URL
    CLINVAR_BASE_URL: str = os.getenv(
        "CLINVAR_BASE_URL", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    )

    # API reliability settings
    GNOMAD_TIMEOUT: int = int(os.getenv("GNOMAD_TIMEOUT", "30"))  # seconds
    GNOMAD_MAX_RETRIES: int = int(os.getenv("GNOMAD_MAX_RETRIES", "3"))
    GNOMAD_RETRY_DELAY: float = float(os.getenv("GNOMAD_RETRY_DELAY", "1.0"))  # seconds

    CLINVAR_MAX_RETRIES: int = int(os.getenv("CLINVAR_MAX_RETRIES", "3"))
    CLINVAR_RETRY_DELAY: float = float(
        os.getenv("CLINVAR_RETRY_DELAY", "1.0")
    )  # seconds

    @classmethod
    def validate(cls) -> None:
        """Validate that required configuration is present."""
        if not cls.NCBI_API_KEY:
            import warnings

            warnings.warn(
                "NCBI_API_KEY not set. API requests will be rate-limited to 3/second. "
                "Register for a free API key at https://www.ncbi.nlm.nih.gov/account/ "
                "and add it to .env file.",
                UserWarning,
            )


# Validate configuration on import
Config.validate()
