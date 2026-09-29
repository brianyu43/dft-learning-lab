"""Join the published DIF pilot observations with explicit 2D descriptors.

This is source curation, not local ML inference or a DFT ranking. Source files
are pinned; an existing differing output must be reviewed as a new version.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import Lipinski, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
SOURCE_HASHES = {
    "data/literature/week1/deepcocrystal.xml": "a7694890b651ee8e49865f46f63a26da52b5360cb21a95d5125b173f1f3f64b3",
    "data/week3/literature/ANIE-64-e202507835-s001.pdf": "684edf33c16e0641481c50d6dcf4eda6af7942f45189f31989d7988c60a2dc0d",
    "data/week3/literature/journal-si-extracted.txt": "44d518afeec56a39ccb661394de8083bfba78525c7f93f1c31b9d83f20f5d510",
    "experiments/week3/literature-review.json": "607ab5cc4ae125e350869f6f426ef6a241158a498fee029a27fa35c0a5fd9b91",
    "reports/week1-molecule-identities.json": "deb13b01ba66b148b8ed1e46e8848f9355e973729de9e9404cfc37a0b958d4ca",
    "experiments/week2/chemistry-candidates.json": "a6bc5fa4ae8bad3488b038a96437aa66a019469f37c217ec4d95833a7b6686b0",
}
MODEL_NAMES = ["DeepCocrystal", "SMI-CNN(can)", "CCGNet", "CC-Descriptor-ML",
               "Fingerprint-DNN", "Descriptor-DNN", "MC", "HBP"]
LABEL = {"✓": "positive", "×": "negative", "?": "abstain", "n.a.": "not_tested"}
KO = {"positive": "형성", "negative": "미형성", "abstain": "보류"}
NAMES = {"adenine": "아데닌", "caffeine": "카페인", "theobromine": "테오브로민", "xanthine": "잔틴"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(condition, message):
    if not condition:
        raise ValueError(message)


def load(rel):
    return json.loads((ROOT / rel).read_text())


def text(node):
    return " ".join("".join(node.itertext()).split())


def published_tables():
    root = ET.parse(ROOT / "data/literature/week1/deepcocrystal.xml")
    wrap = next(w for w in root.iter("table-wrap") if w.get("id") == "anie202507835-tbl-0003")
    headers = [text(c).replace("‐", "-") for c in wrap.findall("./table/thead/tr")[-1]]
    check(headers == ["Coformer", "Experimental", *MODEL_NAMES], "Table 3 header mismatch")
    rows = []
    for row in wrap.findall("./table/tbody/tr"):
        cols = [text(c) for c in row]
        if cols[0].lower() not in {"adenine", "caffeine", "theobromine", "xanthine", "rosmarinic acid", "riboflavin"}:
            continue
        check(len(cols) == 10, "Table 3 row width mismatch")
        rows.append({"coformer": cols[0].lower(), "experimental": LABEL[cols[1]],
                     "published_predictions": dict(zip(MODEL_NAMES, map(LABEL.__getitem__, cols[2:]))),
                     "raw_symbols": cols[1:]})
    check(len(rows) == 6, "Table 3 must retain all six tested pairs")

    si_text = (ROOT / "data/week3/literature/journal-si-extracted.txt").read_text()
    table = si_text.split("PDF PAGE 15\n", 1)[1].split("CO-CRYSTAL CHARACTERIZATION", 1)[0]
    pattern = re.compile(r"^(.+?)\s+(\d\.\d+)\s+(\d\.\d+)\s+±\s+(\d\.\d+)\s+(\d+)%\s+(✓|×|n\.a\.)$")
    si_rows = []
    for line in table.splitlines():
        m = pattern.fullmatch(line.strip())
        if m:
            coformer, smi, mean, sd, vote, outcome = m.groups()
            si_rows.append({"coformer": coformer.lower(), "SMI_CNN_can_score": float(smi),
                            "DeepCocrystal_mean": float(mean), "DeepCocrystal_sd": float(sd),
                            "majority_agreement_percent": int(vote), "experimental": LABEL[outcome],
                            "printed_score_and_sd": f"{mean} ± {sd}", "source_pdf_page": 15})
    check(len(si_rows) == 12 and len({r['coformer'] for r in si_rows}) == 12, "Table S11 requires 12 unique rows")
    by_name = {r['coformer']: r for r in si_rows}
    for row in rows:
        check(row['experimental'] == by_name[row['coformer']]['experimental'], "Table 3 / S11 outcome mismatch")

    # Independent transcription from the visually inspected printed PDF p15.
    reviewed = {"adenine": (.99, .00, 100, "positive"), "caffeine": (.99, .01, 100, "positive"),
                "theobromine": (.66, .35, 60, "negative"), "xanthine": (.63, .38, 60, "negative")}
    for name, expected in reviewed.items():
        row = by_name[name]
        check(tuple(row[k] for k in ['DeepCocrystal_mean', 'DeepCocrystal_sd', 'majority_agreement_percent', 'experimental']) == expected,
              f"Visual transcription mismatch: {name}")

    # The article has different threshold statements. Preserve both; do not
    # reconstruct reported labels by arbitrarily choosing one of them.
    caption = text(wrap.find("caption"))
    body = [text(p) for p in root.iter('p')]
    threshold_paragraphs = [p for p in body if 'For prospective application' in p and 'threshold of 0.10' in p]
    check(len(threshold_paragraphs) == 1 and '0.30' in caption, "Re-review source uncertainty policy")
    return rows, si_rows, {
        "main_text_prospective_sd_threshold": .10,
        "table_3_caption_abstains_when_sd_greater_than": .30,
        "reported_table_labels_are_authoritative_for_this_comparison": True,
        "impact_on_four_pilot_labels": "none; both positive SDs <= 0.10 and both abstained SDs > 0.30",
        "source_locations": ["main text Uncertainty estimation final paragraph", "Table 3 caption"],
        "sd_scope": "variation across augmented SMILES model predictions, not experimental error or a confidence interval",
        "majority_agreement_scope": "agreement with the majority class, not necessarily the positive class and not experimental yield",
        "printed_zero_sd_is_rounded_not_proof_of_zero_variation": True,
    }


def descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    check(mol is not None and len(Chem.GetMolFrags(mol)) == 1, "Invalid or multicomponent descriptor input")
    check(Chem.GetFormalCharge(mol) == 0, "This baseline is explicitly neutral")
    hbd, hba = Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)
    return {"input_smiles": smiles, "canonical_smiles": Chem.MolToSmiles(mol),
            "formula": rdMolDescriptors.CalcMolFormula(mol), "InChIKey": Chem.MolToInchiKey(mol),
            "formal_charge": 0, "HBD_sites": hbd, "HBA_sites": hba,
            "definition_audit": {
                "legacy_python_HDonorSmarts_match_count_not_used": len(mol.GetSubstructMatches(Lipinski.HDonorSmarts)),
                "legacy_python_HAcceptorSmarts_match_count_not_used": len(mol.GetSubstructMatches(Lipinski.HAcceptorSmarts)),
                "atom_indices_not_assigned": "Python SMARTS helpers are not guaranteed equivalent to the compiled CalcNumHBD/HBA descriptors",
            }}


def build():
    for rel, digest in SOURCE_HASHES.items():
        check(sha(ROOT / rel) == digest, f"Source changed; review before regenerating: {rel}")
    check(rdBase.rdkitVersion == "2026.03.6", "Descriptor library changed; version the new baseline")
    all_models, all_scores, uncertainty = published_tables()
    published = {r['coformer']: r for r in all_models}
    scores = {r['coformer']: r for r in all_scores}
    candidates = load("experiments/week2/chemistry-candidates.json")
    states_in_dimer_inventory = {s for r in candidates['dimers'] for s in r['state_ids']}
    features, all_tautomers, identity_source_hashes = {}, [], {}
    for record in load("reports/week1-molecule-identities.json")['molecules']:
        check(record['identity_checks'] == 'passed', "Molecular identity not verified")
        rel = record['source_path']
        check(sha(ROOT / rel) == record['source_sha256'], "PubChem identity source changed")
        identity_source_hashes[rel] = record['source_sha256']
        desc = descriptors(record['pubchem_state_smiles'])
        check(desc['formula'] == record['formula'] and desc['InChIKey'] == record['InChIKey'], "Identity mismatch")
        features[record['name']] = {"molecule": record['name'], "PubChem_CID": record['pubchem_cid'], **desc}
        for state in record['tautomer_inventory']:
            desc_state = descriptors(state['smiles'])
            check(desc_state['formula'] == record['formula'], "Tautomer formula mismatch")
            all_tautomers.append({"molecule": record['name'], "state_id": state['state_id'],
                                 "pubchem_representation": state['is_pubchem_representation'],
                                 "present_in_dimer_candidate_inventory": state['state_id'] in states_in_dimer_inventory,
                                 **desc_state})
    a = features['diflunisal']
    cases = []
    for case in load("experiments/week3/literature-review.json")['cases']:
        name = case['coformer']; b = features[name]; p = published[name]; s = scores[name]
        expected_outcome = 'positive' if case['outcome'] == 'formed' else 'negative'
        check(p['experimental'] == expected_outcome, "Curated experiment / published outcome mismatch")
        check(p['published_predictions']['DeepCocrystal'] == ('positive' if name in ('adenine', 'caffeine') else 'abstain'),
              "Published DeepCocrystal decision changed")
        cases.append({"case_id": case['case_id'], "coformer": name,
                      "literature": case, "table_S11_scores": s,
                      "table_3_published_predictions": p['published_predictions'],
                      "neutral_reference_features": {"A_HBD": a['HBD_sites'], "A_HBA": a['HBA_sites'],
                                                     "B_HBD": b['HBD_sites'], "B_HBA": b['HBA_sites'],
                                                     "at_least_one_cross_donor_acceptor_direction": bool(
                                                         a['HBD_sites'] * b['HBA_sites'] + b['HBD_sites'] * a['HBA_sites'])},
                      "calculation_condition_comparison": {
                          "baseline_environment": "isolated vacuum molecule/pair",
                          "planned_sensitivity_environment": "SMD ethanol",
                          "ethanol_documented_for_this_pair_preparation": True if 'preparation' in case else None,
                          "full_experimental_condition_match": False,
                          "negative_attempt_conditions": case.get('attempt_conditions'),
                          "local_pair_composition": "1 molecule A + 1 molecule B",
                          "crystallographic_product_stoichiometry_verified": False,
                      },
                      "DFT_B_kcal_mol": None, "DFT_X_kcal_mol": None,
                      "DFT_score_status": "not_assembled; common numerical protocol and full control panel pending"})
    check(len(cases) == 4, "Pilot denominator must remain four pairs")
    check(all(c['neutral_reference_features']['at_least_one_cross_donor_acceptor_direction'] for c in cases),
          "Re-review the reported descriptor observation")
    return {
        "artifact_version": "dif-pilot-baselines-r1", "review_date": "2026-09-29",
        "status": "source_and_descriptor_baseline_complete_DFT_comparison_pending",
        "source_doi": "10.1002/anie.202507835", "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12258686/",
        "source_sha256": {**SOURCE_HASHES, **identity_source_hashes},
        "visual_review": {"pdf_pages": [3, 15], "table": "S11", "outcomes_and_conditions_checked": True},
        "script_sha256": sha(Path(__file__)),
        "descriptor_method": {"rdkit_version": rdBase.rdkitVersion,
            "functions": ["rdkit.Chem.Lipinski.NumHDonors", "rdkit.Chem.Lipinski.NumHAcceptors"],
            "donor_definition_version": Lipinski.NumHDonors.version,
            "acceptor_definition_version": Lipinski.NumHAcceptors.version,
            "Lipinski_python_sha256": sha(Path(Lipinski.__file__)),
            "rdMolDescriptors_binary_sha256": sha(Path(rdMolDescriptors.__file__)),
            "primary_representation": "validated neutral PubChem graph; not an optimized or proton-transferred complex",
            "counts_are": "donor/acceptor sites in a 2D graph, not donor hydrogens or realized hydrogen bonds",
            "not_included": ["3D accessibility", "intramolecular hydrogen bonding", "bond strength", "crystal packing", "proton-transfer populations"],
            "is_HBP_model": False, "class_predictor_fitted": False},
        "uncertainty_policy_source_discrepancy": uncertainty,
        "published_predictions_origin": "paper tables, not local model execution",
        "published_model_evaluation_design": {"DeepCocrystal": "prospective in original study; retrospective reuse here",
                                               "other_Table_3_models": "retrospective in original study"},
        "paper_table_3_all_six_tested_pairs": all_models,
        "paper_table_S11_all_twelve_candidates": all_scores,
        "neutral_reference_molecule_features": list(features.values()),
        "enumerated_tautomer_features": all_tautomers,
        "tautomer_inventory_scope": "enumerated graphs, not measured populations, stable states, or polymorphs",
        "pilot_cases": cases,
        "pilot_denominators": {"molecular_pairs": 4, "observed_formed": 2, "observed_not_formed": 2,
                               "published_DeepCocrystal_decisions": 2, "published_DeepCocrystal_abstentions": 2,
                               "independent_experimental_replicates": None},
        "limitations": [
            "Known-outcome selected pilot; not an independent model performance test.",
            "Missing negative attempt conditions remain unknown, not assumed vacuum, ethanol, room temperature, or 1:1.",
            "Published spectral proton-transfer interpretation is not a solved crystal structure or polymorph inventory.",
            "The existence of a cross donor/acceptor direction is shared by all four pairs; no fitted threshold or new accuracy claim.",
            "DFT and local ML results are absent; no zero-value substitution for missing scores.",
        ],
    }


def render_report(d):
    lines = ["# 디플루니살 네 사례의 문헌·모델·수소결합 특징 비교", "",
        "2026-09-29. 정식 논문 Table 3, 정식 보충자료 Table S11 및 실험법을 대조한 기준표다. "
        "DFT 점수와 로컬 모델 추론은 아직 포함하지 않았다. 현재 수치 정밀도와 전체 대조군 검증을 마친 뒤 같은 사례에 B/X를 추가한다.", "",
        "## 실험 조건과 계산 환경", "",
        "| 사례 | 문헌 관찰 | 실제 보고된 제조 조건 | 아직 모르는 것 |", "|---|---|---|---|",
        "| D01 DIF–아데닌 | 다성분 결정 형성 | DIF 100 mg + ADE 54 mg, 투입 1:1, 에탄올 1 mL, 상온 슬러리 48시간 후 용매 증발 | 정확한 온도·반복 수·회절로 확정한 구조와 화학량론 |",
        "| D02 DIF–카페인 | 다성분 결정 형성 | DIF 100 mg + CAF 78 mg, 투입 1:1, 건식 분쇄 20분 + 에탄올 보조 분쇄 20분; 5분마다 에탄올 추가 | 에탄올 총량·온도·반복 수·회절로 확정한 구조와 화학량론 |",
        "| D03 DIF–테오브로민 | 해당 실험에서 미형성 | 시도별 조건 미보고 | 용매·몰비·온도·시간·반복 수 등 |",
        "| D04 DIF–잔틴 | 해당 실험에서 미형성 | 시도별 조건 미보고 | 용매·몰비·온도·시간·반복 수 등 |", "",
        "형성 두 사례는 FT-IR ATR 및 고체 13C/15N NMR로 분석했다. 문헌은 CAF에서 COOH···N9, "
        "ADE에서 COO−···+HN1 및 N10이 참여하는 접촉을 해석한다. ADE에는 양성자 이동 가능성이 있어 중성 분자의 특징 개수와 실제 이온쌍을 구별한다. "
        "이 해석은 원자 좌표나 결정 다형 수를 제공하지 않는다.", "",
        "우리의 기본 계산은 진공의 고립 분자/분자쌍이며, 에탄올 SMD는 환경 민감도 확인이다. "
        "고체 분쇄·슬러리·증발·결정 패킹을 재현한 계산이 아니다. 특히 미형성 두 사례의 용매를 에탄올로 채우지 않는다. "
        "SI의 일반 용매 목록과 외부 데이터셋의 반복 횟수도 이 네 사례의 개별 시도에 복사하지 않는다.", "",
        "## 논문에 보고된 예측", "",
        "| 상대 분자 | 실험 | DeepCocrystal 점수 평균 ± SD | 보고된 최종 판정 | 다수결 일치율 | SMI-CNN(can) 점수 |",
        "|---|---|---:|---|---:|---:|"]
    for c in d['pilot_cases']:
        name=c['coformer']; s=c['table_S11_scores']; p=c['table_3_published_predictions']
        lines.append(f"| {NAMES[name]} | {KO[s['experimental']]} | {s['printed_score_and_sd']} | {KO[p['DeepCocrystal']]} | {s['majority_agreement_percent']}% | {s['SMI_CNN_can_score']:.2f} |")
    lines += ["", "SD는 증강 SMILES에 따른 예측 변동이다. 실험의 오차막대·반복 성공률·신뢰구간이 아니다. "
        "0.00은 표시 자릿수에서 반올림된 값이다. 테오브로민·잔틴의 평균이 0.5보다 높아도, 논문의 최종 판정은 보류다. "
        "우리의 네 쌍 중 판정 2쌍, 보류 2쌍이다. 논문 전체 여섯 쌍의 성능 수치를 이 부분집합으로 가져오지 않는다.", "",
        "**출처 내 불일치:** 본문은 전향 적용의 SD 기준을 0.10으로, Table 3 설명은 보류 기준을 SD > 0.30으로 적는다. "
        "두 문구를 모두 보존하고 보고된 최종 판정을 사용했다. 이 네 사례에서는 양성 두 SD가 0.10 이하, 보류 두 SD가 0.30 초과여서 판정 차이는 없다.", "",
        "| 상대 분자 | CCGNet | CC-Descriptor-ML | Fingerprint-DNN | Descriptor-DNN | MC | HBP |",
        "|---|---|---|---|---|---|---|"]
    for c in d['pilot_cases']:
        p=c['table_3_published_predictions']
        lines.append('| '+NAMES[c['coformer']]+' | '+' | '.join(KO[p[m]] for m in MODEL_NAMES[2:])+' |')
    lines += ["", "이 값은 저자들이 Table 3에 보고한 결과다. 여기서 해당 모델을 실행하지 않았다. "
        "원논문에서 DeepCocrystal은 전향 예측이고 다른 모델들은 사후 비교이며, 우리 프로젝트는 알려진 결과를 재사용하는 사후 파일럿이다.", "",
        "## 직접 계산한 단순 수소결합 특징", "",
        "RDKit 2026.03.6의 `Lipinski.NumHDonors` / `NumHAcceptors`를 검증된 중성 PubChem 표현에 적용했다. "
        "HBD/HBA는 주는 자리/받는 자리의 개수이며, 실제 형성된 수소결합 수나 세기가 아니다. "
        "디플루니살 자체는 HBD 2, HBA 2다. 구현 정의가 다른 단순 N/O 개수나 초기 자세 생성기의 접촉 후보 수와 섞지 않는다.", "",
        "| 상대 분자 | HBD | HBA | DIF와 주는 자리–받는 자리 조합 존재 | 문헌 관찰 |",
        "|---|---:|---:|---|---|"]
    for c in d['pilot_cases']:
        f=c['neutral_reference_features']
        lines.append(f"| {NAMES[c['coformer']]} | {f['B_HBD']} | {f['B_HBA']} | 있음 | {KO[c['table_S11_scores']['experimental']]} |")
    lines += ["", "네 쌍 모두 상대 수소결합 자리를 갖는다. 따라서 '그런 자리가 존재한다'는 특징 하나는 이 네 사례를 구분하지 못한다. "
        "개수 전체가 무용하다는 뜻은 아니며, 알려진 네 라벨에 맞춰 개수 임계값을 학습하지 않았다. "
        "3차원 접근성, 분자 내부 수소결합, 양성자 이동, 결정 패킹은 이 특징에 포함되지 않는다. "
        "이 계산은 표 위의 CCDC Hydrogen Bond Propensity(HBP) 모델과 다르다.", "",
        "구현 점검에서 RDKit의 별도 Python `HAcceptorSmarts` 패턴은 CAF/THEO를 각각 6/5개로 세지만, "
        "지정한 `NumHAcceptors`가 호출하는 `CalcNumHBA`는 각각 3/3개로 센다는 차이를 확인했다. "
        "주 표는 명시한 함수의 출력으로 통일하고, 두 정의를 혼합한 원자별 대응은 만들지 않았다. 차이는 JSON에 별도로 보존한다.", "",
        f"JSON에는 원래 열거한 {len(d['enumerated_tautomer_features'])}개 토토머 표현의 특징도 보관했다. "
        "이는 안정한 상태나 실제 존재 비율, 다형 개수의 측정이 아니다. 주 비교의 중성 표현과 후보 목록에 포함된 대안 표현을 구별한다.", "",
        "## 다음 DFT 비교와 남은 요구사항", "",
        "B = E(AB) − E(A) − E(B), X = E(AB) − ½E(AA) − ½E(BB)를 동일한 방법·환경에서 검증된 최소점으로 구성한다. "
        "현재 두 점수 열은 JSON/CSV에서 미산출로 남긴다. 단순한 접촉 가능성에서 더 나아가 자기 결합과의 경쟁을 고려했을 때, "
        "형성 2쌍과 미형성 2쌍의 교차 비교 네 건을 더 일관되게 설명하는지 평가한다. 동률 범위는 기존 계획의 0/0.5/1/2 kcal/mol을 모두 사용한다.", "",
        "D06의 문헌/기존 모델/2D 특징 부분은 준비됐지만 DFT와의 실제 비교는 미완료다. "
        "개별 미형성 조건을 확보하지 못한 한 그 실험이 실패한 원인을 특정하지 않는다. "
        "새 독립 사례로 일반화하거나 기존 모델보다 낫다고 결론 내리지 않는다.", "",
        "## 출처와 재현", "",
        "- [정식 논문 Table 3 및 본문](https://pmc.ncbi.nlm.nih.gov/articles/PMC12258686/)",
        "- [정식 보충자료](../../data/week3/literature/ANIE-64-e202507835-s001.pdf): PDF 3쪽 실험법, 15쪽 Table S11 및 해석; 두 페이지 렌더링 재검토.",
        "- [세부 JSON](drug-pilot-baselines-r1.json), [네 사례 CSV](drug-pilot-baselines-r1.csv): 원문 12후보/6시험 전체, 4사례, 구조 표현, 출처 해시와 소프트웨어 정의 포함.",
        "- 재생성: `.venv/bin/python scripts/curate_completion_drug_baselines.py`. 기존 파일이 다르면 덮어쓰지 않는다.", ""]
    return '\n'.join(lines)


def main():
    data = build()
    out = ROOT/'reports/completion/drug-pilot-baselines-r1.json'
    rows = []
    for c in data['pilot_cases']:
        s=c['table_S11_scores']
        rows.append({"case_id":c['case_id'], "coformer":c['coformer'], "experimental":s['experimental'],
            "DeepCocrystal_mean":s['DeepCocrystal_mean'], "DeepCocrystal_sd":s['DeepCocrystal_sd'],
            **c['table_3_published_predictions'], **c['neutral_reference_features'],
            "DFT_B_kcal_mol":None, "DFT_X_kcal_mol":None, "DFT_score_status":c['DFT_score_status']})
    buf=io.StringIO(); writer=csv.DictWriter(buf,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    outputs = {out: json.dumps(data,ensure_ascii=False,indent=2)+'\n', out.with_suffix('.csv'):buf.getvalue(),
               ROOT/'reports/completion/DRUG_PILOT_BASELINES.ko.md':render_report(data)}
    for p, content in outputs.items():
        check(not p.exists() or p.read_bytes() == content.encode(), f"Output exists with different contents: {p}")
    for p, content in outputs.items():
        p.parent.mkdir(parents=True,exist_ok=True)
        if not p.exists():p.write_bytes(content.encode())
    print(json.dumps({"status":data['status'], "cases":len(rows), "tautomer_graphs":len(data['enumerated_tautomer_features']),
                      "outputs":{str(p.relative_to(ROOT)):sha(p) for p in outputs}},indent=2))


if __name__ == '__main__':
    main()
