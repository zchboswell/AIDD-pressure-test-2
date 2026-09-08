# CEDAR input package

Read `REQUEST.md`, then the project material you need. All inputs arrive together; there are no later rounds or surprise messages. Earlier and current records coexist as they would in a project handover.

The `data/` directory contains the activity archive, inventory, simulated program readouts, vendor features, and experiment menu. The `references/` directory contains selected public source records and assay metadata. The `structures/` directory contains the deposited experimental object and derived vendor exports. The `legacy/` directory contains a small analysis helper and notebook-style model script. The `project/` directory contains scope, provenance, and inherited experience.

The original files are read-only inputs. Implement repairs and derived data in your repository. Review supplied scripts before running them. No serialized third-party model is supplied. Use a small reproducible model or another justified method if useful; there is no requirement to train a complex model or to dock compounds.

Source records are a selected snapshot, not a complete ChEMBL download. Some prospective-in-this-simulation candidate measurements are intentionally absent. Do not retrieve candidate activities externally. Public chemical identifiers are retained for source fidelity; no claim is made that these molecules are unseen in model pretraining.
