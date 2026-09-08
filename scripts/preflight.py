import hashlib, importlib, json, pathlib, datetime
root = pathlib.Path(__file__).resolve().parents[1]
manifest_path = root / 'inputs/MANIFEST.json'
manifest = json.loads(manifest_path.read_bytes())
checks = []
for name, expected in manifest['files'].items():
    path = root / 'inputs' / name
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    checks.append({'file': name, 'sha256': actual, 'matches': actual == expected})
env = json.loads((root / 'admin/ENVIRONMENT.json').read_text())
packages = {}
for name, expected in env['packages'].items():
    try:
        mod = importlib.import_module(name)
        version = str(mod.cmd.get_version()[0]) if name == 'pymol' else str(mod.__version__)
        packages[name] = {'import_ok': True, 'version': version, 'expected_version': expected['version'], 'version_matches': version == expected['version']}
    except Exception as exc:
        packages[name] = {'import_ok': False, 'error': str(exc)}
skills = ['imagegen','openai-docs','plugin-creator','skill-creator','skill-installer','method-capability-review','potency-modeling','scientific-software','structural-analysis','breakdown','deep-research-work:deep-research','documents:documents','pdf:pdf','plugin-management:plugin-management','presentations:Presentations','qwen-worker:qwen-ops-review','sites:sites-building','sites:sites-hosting','spreadsheets:Spreadsheets','spreadsheets:excel-live-control','template-creator:template-creator','visualize:visualize']
record = {'phase':'operational_preflight', 'checked_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(), 'assigned_model':'GPT-6 Astra', 'assigned_effort':'Medium', 'observed_model':None, 'observed_effort':None, 'observation_note':'Assignment is explicit; no independent runtime model/effort telemetry available.', 'offered_skill_names':skills, 'guidance_consistency':'Four offered scientific skills match the assigned AIDD-equipped arm. Generic deep-research/internet capabilities remain constrained by the closed scientific internet boundary.', 'substantive_skill_files_read':[], 'scientific_input_content_read':False, 'manifest_sha256':hashlib.sha256(manifest_path.read_bytes()).hexdigest(), 'file_checks':checks, 'all_manifest_hashes_match':all(c['matches'] for c in checks), 'runtime_command':env['command'], 'runtime_checks':packages, 'private_repo_url':None, 'private_visibility_verified':False, 'commit':None, 'operational_blocker':None, 'scientific_work_started':False}
(root / 'admin/PREFLIGHT.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({'manifest_sha256':record['manifest_sha256'],'files_verified':len(checks),'all_match':record['all_manifest_hashes_match'],'runtime_checks':packages},indent=2))
