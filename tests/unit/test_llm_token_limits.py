from utils.llm_utils import clamp_max_tokens


def test_new_azure_deployments_have_large_output_budget():
    assert clamp_max_tokens("azure_ai/gpt-5.4", 60000, warn=False) == 15000
    assert clamp_max_tokens("azure_ai/DeepSeek-V4-Pro", 60000, warn=False) == 15000
