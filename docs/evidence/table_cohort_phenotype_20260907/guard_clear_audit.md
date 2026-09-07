# Guard-clear audit on the two locked cont120 candidates (2026-09-07)

What the always-on phenotype guard cleared, restored into a copy of the locked predictions (fill-null-only) and rescored. Zero LLM calls. Verdict 'exact' means the raw cleared value equals gold; 'wrong' means restoring it creates a new count error.

```
== guard clears restored: by (field, reason, verdict) ==
   10 ('affected', 'copied_carriers_onto_affected', 'exact')
    6 ('affected', 'copied_carriers_onto_affected', 'wrong')
    1 ('affected', 'copied_carriers_onto_affected', 'unmatched')
   77 ('affected', 'unsourced_zero_affected', 'exact')
    1 ('affected', 'unsourced_zero_affected', 'unmatched')
   27 ('unaffected', 'implied_unaffected_zero', 'exact')
    2 ('unaffected', 'implied_unaffected_zero', 'wrong')

== examples ==
('affected', 'unsourced_zero_affected', 'exact')
     ('KCNH2', '29331839', 'p.Pro963Thr', 0, (2, 0, 2), 'llm_text', 'Results, Clinical and genetic ')
     ('KCNQ1', '17038145', 'p.Gly589Asp', 0, (9, 0, 9), 'llm_text', 'Methods, Patient Population; R')
     ('SCN5A', '20129283', 'A286S c.856G>T', 0, (1, 0, 1), 'regex_table', 'Table 2')
     ('SCN5A', '20129283', 'A447G c.1340C>G', 0, (1, 0, 1), 'regex_table', 'Table 2')
('affected', 'copied_carriers_onto_affected', 'exact')
     ('RYR2', '32218223', 'p.Asp2216Gly', 2, (2, 2, 0), 'llm_text', 'Results, section 3.3 Genetic S')
     ('SCN5A', '12051963', 'p.Gly351Val', 3, (3, 3, 0), 'llm_text', 'Results, Mutation analysis, Fa')
     ('KCNQ1', '30036649', 'p.Ala344= c.1032G>C', 4, (6, 4, 0), 'llm_table', 'Table 2, Results section Genet')
     ('KCNQ1', '30036649', 'p.Gly314Ser c.940G>A', 2, (2, 2, 0), 'llm_table', 'Table 2, Results Genetic analy')
('unaffected', 'implied_unaffected_zero', 'exact')
     ('RYR2', '32218223', 'p.Asp2216Gly', 0, (2, 2, 0), 'llm_text', 'Results, section 3.3 Genetic S')
     ('KCNH2', '21130771', 'p.T618I c.1853C>T', 0, (0, 0, 0), 'llm_text', 'Results, mutation analysis; Ab')
     ('KCNQ1', '30036649', 'p.Ala344= c.1032G>C', 0, (6, 4, 0), 'llm_table', 'Table 2, Results section Genet')
     ('KCNQ1', '30036649', 'p.Gly314Ser c.940G>A', 0, (2, 2, 0), 'llm_table', 'Table 2, Results Genetic analy')
('affected', 'copied_carriers_onto_affected', 'wrong')
     ('KCNH2', '21130771', 'p.T618I c.1853C>T', 4, (0, 0, 0), 'llm_text', 'Results, mutation analysis; Ab')
     ('KCNH2', '17171344', 'p.Gly604Ser c.1810G>A', 11, (10, 9, 1), 'llm_text', 'Abstract')
     ('SCN5A', '15338453', 'p.Gly1262Ser c.3934G>A', 4, (2, 2, 0), 'llm_text', 'Abstract')
     ('SCN5A', '28584071', 'R1193Q c.3578G>A', 6, (0, 0, 0), 'llm_table', 'Results, Tables 1 and 3')
('unaffected', 'implied_unaffected_zero', 'wrong')
     ('KCNH2', '17171344', 'p.Gly604Ser c.1810G>A', 0, (10, 9, 1), 'llm_text', 'Abstract')
     ('KCNQ1', '28491547', 'p.Val141Met c.421G>A', 0, (6, 3, 3), 'llm_text', 'Results, Genetic analysis')
('affected', 'unsourced_zero_affected', 'unmatched')
     ('KCNQ1', '19841300', 'F335L', 0, None, 'regex_table', 'Table 2. Compendium Summary of')
('affected', 'copied_carriers_onto_affected', 'unmatched')
     ('SCN5A', '28584071', 'p.His558Arg c.1673A>G', 9, None, 'llm_table', 'Results, Table 1')

== papers with >=3 restored ==
('SCN5A', '20129283') {('affected', 'exact'): 55}
('KCNQ1', '19841300') {('affected', 'exact'): 20, ('unaffected', 'exact'): 18, ('affected', 'unmatched'): 1}
('KCNQ1', '30036649') {('affected', 'exact'): 2, ('unaffected', 'exact'): 2}
('RYR2', '12093772') {('affected', 'wrong'): 2, ('unaffected', 'exact'): 2}
```

Reading: relaxing the copy guard for model rows would recover 10 exact affected values and add 6 wrong ones; the zero restorations land only on zero-gold rows (no positive-gold recovery). Not a lever; the table-cohort projection recovers the deterministic-table cases without touching model rows.
