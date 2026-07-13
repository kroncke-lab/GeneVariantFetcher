from harvesting.elsevier_api import ElsevierAPIClient
from pipeline.extraction import ExpertExtractor


def test_elsevier_doi_detection_excludes_nature_prefix():
    assert ElsevierAPIClient.is_elsevier_doi("10.1016/j.hrthm.2018.01.014") is True
    assert ElsevierAPIClient.is_elsevier_doi("10.1038/s41467-021-27599-6") is False


def test_elsevier_xml_to_markdown_preserves_table_bodies_for_deterministic_parser():
    xml = """
    <full-text-retrieval-response>
      <coredata><dc:title>Compendium</dc:title></coredata>
      <originalText>
        <body>
          <section>
            <section-title>Results</section-title>
            <para>Overall, KCNH2 and SCN5A variants were identified.</para>
          </section>
        </body>
      </originalText>
      <table id="tbl3">
        <label>Table 3</label>
        <caption><simple-para>Summary of KCNH2 LQT2-associated variants</simple-para></caption>
        <tgroup cols="6">
          <thead>
            <row><entry>No.</entry><entry>Exon</entry><entry>Nucleotide</entry><entry>Variant</entry><entry>Location</entry><entry>No. of patients</entry></row>
          </thead>
          <tbody>
            <row><entry>27</entry><entry>6</entry><entry>1264 G&#x2192;A</entry><entry>A422T &#x204E;</entry><entry>S1</entry><entry>1</entry></row>
            <row><entry>28</entry><entry>6</entry><entry>del C 287</entry><entry>S95fs/141 &#x204E;</entry><entry>N-terminus</entry><entry>1</entry></row>
            <row><entry>51</entry><entry>7</entry><entry>1904 A&#x2192;T</entry><entry>N6351 &#x204E;</entry><entry>PORE-S6</entry><entry>1</entry></row>
            <row><entry>61</entry><entry>9</entry><entry>2398+5 G&#x2192;T</entry><entry>L799/sp &#x204E;</entry><entry>C-terminus</entry><entry>3</entry></row>
            <row><entry>73</entry><entry>12</entry><entry>2143+5 G&#x2192;A</entry><entry>A715sp</entry><entry>C-terminus</entry><entry>1</entry></row>
            <row><entry>68</entry><entry>10</entry><entry>2592+3 G&#x2192;A</entry><entry>D864/sp &#x204E;</entry><entry>C-terminus</entry><entry>1</entry></row>
          </tbody>
        </tgroup>
      </table>
      <table id="tbl4">
        <label>Table 4</label>
        <caption><simple-para>Summary of SCN5A LQT3-associated variants</simple-para></caption>
        <tgroup cols="6">
          <thead>
            <row><entry>No.</entry><entry>Exon</entry><entry>Nucleotide</entry><entry>Variant</entry><entry>Location</entry><entry>No. of patients</entry></row>
          </thead>
          <tbody>
            <row><entry>20</entry><entry>26</entry><entry>del 4511-4519</entry><entry>KPQ1505-1507del</entry><entry>IIIS6-IVS1</entry><entry>1</entry></row>
            <row><entry>28</entry><entry>28</entry><entry>5350 G&#x2192;A</entry><entry>E1784K</entry><entry>C-terminus</entry><entry>4</entry></row>
          </tbody>
        </tgroup>
      </table>
    </full-text-retrieval-response>
    """

    markdown = ElsevierAPIClient(api_key="key").xml_to_markdown(xml)

    assert markdown is not None
    assert "Table 3 Summary of KCNH2 LQT2-associated variants" in markdown
    assert "| 27 | 6 | 1264 G>A | A422T ⁎ | S1 | 1 |" in markdown

    extractor = ExpertExtractor(models=["gpt-4"])
    kcnh2_variants = extractor._parse_markdown_table_variants(markdown, "KCNH2")
    by_protein = {v["protein_notation"]: v for v in kcnh2_variants}

    assert set(by_protein) == {
        "A422T",
        "S95fsX",
        "N635I",
        "L799sp",
        "A715sp",
        "D864sp",
    }
    assert by_protein["A422T"]["cdna_notation"] == "c.1264G>A"
    assert by_protein["A422T"]["penetrance_data"] == {
        "total_carriers_observed": 1,
        "affected_count": 1,
        "unaffected_count": 0,
    }

    scn5a_variants = extractor._parse_markdown_table_variants(markdown, "SCN5A")
    assert [v["protein_notation"] for v in scn5a_variants] == [
        "K1505_Q1507del",
        "E1784K",
    ]
    assert scn5a_variants[1]["penetrance_data"]["total_carriers_observed"] == 4


def test_extract_supplement_refs_finds_mmc_and_dedups():
    from harvesting.elsevier_api import ElsevierAPIClient

    xml = (
        "<xocs:pii-unformatted>S1547527118300146</xocs:pii-unformatted>"
        '<foo href="https://ars.els-cdn.com/content/image/'
        '1-s2.0-S1547527118300146-mmc1.docx"/>'
        "see also 1-s2.0-S1547527118300146-mmc1.docx and "
        "1-s2.0-S1547527118300146-mmc2.xlsx"
    )
    refs = ElsevierAPIClient.extract_supplement_refs(xml)
    assert refs == [
        "1-s2.0-S1547527118300146-mmc1.docx",
        "1-s2.0-S1547527118300146-mmc2.xlsx",
    ]  # de-duplicated, order-preserving
    assert ElsevierAPIClient.extract_supplement_refs("") == []
    assert ElsevierAPIClient.extract_supplement_refs("no supplements here") == []
