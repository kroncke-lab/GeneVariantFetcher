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


def test_llm_provider_key_status_accepts_anthropic(monkeypatch):
    for key in bootstrap.LLM_PROVIDER_KEY_ENV_VARS:
        monkeypatch.delenv(key, raising=False)

    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-ant")

    assert bootstrap.has_llm_provider_key()
    assert bootstrap.llm_provider_key_status()["ANTHROPIC_API_KEY"] is True


def test_network_hardening_noop_by_default(monkeypatch):
    import socket

    monkeypatch.delenv("GVF_FORCE_IPV4", raising=False)
    monkeypatch.delenv("GVF_SOCKET_TIMEOUT", raising=False)

    getaddrinfo_before = socket.getaddrinfo
    timeout_before = socket.getdefaulttimeout()
    try:
        bootstrap._apply_network_hardening()
        assert socket.getaddrinfo is getaddrinfo_before
        assert socket.getdefaulttimeout() == timeout_before
    finally:
        socket.getaddrinfo = getaddrinfo_before
        socket.setdefaulttimeout(timeout_before)


def test_force_ipv4_resolves_af_inet_only(monkeypatch):
    import socket

    monkeypatch.setenv("GVF_FORCE_IPV4", "1")
    monkeypatch.delenv("GVF_SOCKET_TIMEOUT", raising=False)

    seen_families = []

    def spy(host, port, family=0, type=0, proto=0, flags=0):
        seen_families.append(family)
        return [(socket.AF_INET, socket.SOCK_STREAM, 6, "", (host, port))]

    getaddrinfo_before = socket.getaddrinfo
    # monkeypatch records `getaddrinfo_before` and restores it on teardown, even
    # though _apply_network_hardening reassigns socket.getaddrinfo afterward.
    monkeypatch.setattr(socket, "getaddrinfo", spy)
    try:
        bootstrap._apply_network_hardening()
        socket.getaddrinfo("example.com", 443, socket.AF_INET6)
        assert seen_families[-1] == socket.AF_INET
    finally:
        socket.getaddrinfo = getaddrinfo_before


def test_socket_timeout_sets_process_default(monkeypatch):
    import socket

    monkeypatch.delenv("GVF_FORCE_IPV4", raising=False)
    monkeypatch.setenv("GVF_SOCKET_TIMEOUT", "42")

    timeout_before = socket.getdefaulttimeout()
    try:
        bootstrap._apply_network_hardening()
        assert socket.getdefaulttimeout() == 42.0
    finally:
        socket.setdefaulttimeout(timeout_before)


def test_socket_timeout_ignores_invalid_value(monkeypatch):
    import socket

    monkeypatch.delenv("GVF_FORCE_IPV4", raising=False)
    monkeypatch.setenv("GVF_SOCKET_TIMEOUT", "not-a-number")

    timeout_before = socket.getdefaulttimeout()
    try:
        bootstrap._apply_network_hardening()
        assert socket.getdefaulttimeout() == timeout_before
    finally:
        socket.setdefaulttimeout(timeout_before)
