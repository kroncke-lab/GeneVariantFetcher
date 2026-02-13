"""Tests for runtime bootstrap initialization."""

from utils import bootstrap


def test_initialize_runtime_loads_extra_env_file(monkeypatch, tmp_path):
    env_file = tmp_path / "extra.env"
    env_file.write_text("GVF_BOOTSTRAP_TEST_KEY=loaded\n", encoding="utf-8")

    monkeypatch.setenv("GVF_EXTRA_ENV_FILE", str(env_file))
    monkeypatch.delenv("GVF_BOOTSTRAP_TEST_KEY", raising=False)
    monkeypatch.setattr(bootstrap, "_BOOTSTRAPPED", False)

    bootstrap.initialize_runtime()

    assert bootstrap.os.getenv("GVF_BOOTSTRAP_TEST_KEY") == "loaded"


def test_initialize_runtime_is_idempotent(monkeypatch, tmp_path):
    env_file = tmp_path / "idempotent.env"
    env_file.write_text("GVF_BOOTSTRAP_IDEMPOTENT=first\n", encoding="utf-8")

    monkeypatch.setenv("GVF_EXTRA_ENV_FILE", str(env_file))
    monkeypatch.delenv("GVF_BOOTSTRAP_IDEMPOTENT", raising=False)
    monkeypatch.setattr(bootstrap, "_BOOTSTRAPPED", False)

    bootstrap.initialize_runtime()
    env_file.write_text("GVF_BOOTSTRAP_IDEMPOTENT=second\n", encoding="utf-8")
    bootstrap.initialize_runtime()

    assert bootstrap.os.getenv("GVF_BOOTSTRAP_IDEMPOTENT") == "first"
