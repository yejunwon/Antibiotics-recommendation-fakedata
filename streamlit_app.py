import streamlit as st
import json
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

# ---- 데이터 로딩 ----
with open("Patient Example.json", "r", encoding="utf-8") as f:
    patients = json.load(f)

# 업데이트된 항생제 목록
abx_nodes = [
    "Tazoferan(R) 4.5g",
    "cefaZOLin 1g",
    "Azithromycin 250mg",
    "cefTRIAXone sod 2g",
    "cefePIMe 1g",
    "Amoxclan duo(R) 437.5mg/62.5mg",
    "Meropenem 500mg"
]

# 그람 양성/음성 적용 범위
abx_to_gram = {
    "Tazoferan(R) 4.5g": ["gram_positive", "gram_negative"],  # 광범위
    "cefaZOLin 1g": ["gram_positive"],                        # 주로 MSSA, 일부 GN
    "Azithromycin 250mg": ["gram_positive"],                  # 주로 GP, 비정형균
    "cefTRIAXone sod 2g": ["gram_negative"],                  # 일부 GP 커버 가능하나 주 대상은 GN
    "cefePIMe 1g": ["gram_negative"],                         # GN 위주, P.aeruginosa 포함
    "Amoxclan duo(R) 437.5mg/62.5mg": ["gram_positive", "gram_negative"],  # 광범위, E.faecalis 포함
    "Meropenem 500mg": ["gram_positive", "gram_negative"]     # 광범위, ESBL 포함
}

# 항생제별 독성 가중치(배수). 1.0 = 기준, >1.0 = 더 위험, <1.0 = 상대적으로 안전
# (임상 상식 기반의 예시값. 니 환경에 맞게 쉽게 조정 가능)
abx_risk = {
    "Tazoferan(R) 4.5g": {"age": 1.0, "renal": 1.1, "hepatic": 1.0},
    "cefaZOLin 1g": {"age": 0.9, "renal": 1.0, "hepatic": 0.9},
    "Azithromycin 250mg": {"age": 0.9, "renal": 0.8, "hepatic": 1.2},   # 간대사 고려
    "cefTRIAXone sod 2g": {"age": 1.0, "renal": 1.0, "hepatic": 1.1},   # 담즙배설 이슈 고려
    "cefePIMe 1g": {"age": 1.0, "renal": 1.3, "hepatic": 1.0},          # 신독성 리스크 상대가중
    "Amoxclan duo(R) 437.5mg/62.5mg": {"age": 1.0, "renal": 1.0, "hepatic": 1.1},
    "Meropenem 500mg": {"age": 1.0, "renal": 1.2, "hepatic": 1.0}
}

# === (중요) 0~12 스케일의 환자 기반 baseline toxicity 산출 ===
# 필요 시 너희 병원 스케일로 컷 수정하면 됨.
def baseline_age_toxicity(age: int) -> int:
    # 0~12로 클램프
    if age <= 30:   return 1
    if age <= 50:   return 4
    if age <= 65:   return 7
    if age <= 75:   return 9
    if age <= 85:   return 11
    return 12

def baseline_renal_toxicity(creat: float) -> int:
    # Cr만 있는 예시. eGFR 있으면 그걸로 바꾸는 걸 권장.
    if creat <= 1.0:    return 2
    if creat <= 1.5:    return 4
    if creat <= 2.0:    return 7
    if creat <= 3.0:    return 9
    return 12

def baseline_hepatic_toxicity(bili: float, ast: float=None, alt: float=None) -> int:
    # 단순 Bilirubin 중심 예시. AST/ALT 유무에 따라 강화 가능.
    score = 0
    if bili <= 1.2:     score = 2
    elif bili <= 2.0:   score = 5
    elif bili <= 3.0:   score = 8
    else:               score = 12
    # AST/ALT가 많이 높으면 약간 가산 (선택)
    if ast is not None and alt is not None:
        if ast > 100 or alt > 100:
            score = min(12, score + 2)
    return score

def patient_baseline_toxicity(patient):
    age = patient.get('age', 0)
    creat = patient.get('renal_function', {}).get('creatinine', 0.0)
    bili = patient.get('hepatic_function', {}).get('bilirubin', 0.0)
    ast = patient.get('hepatic_function', {}).get('AST', None)
    alt = patient.get('hepatic_function', {}).get('ALT', None)
    return {
        "age": baseline_age_toxicity(age),
        "renal": baseline_renal_toxicity(creat),
        "hepatic": baseline_hepatic_toxicity(bili, ast, alt)
    }

def abx_specific_toxicity(abx: str, base: dict) -> dict:
    """항생제별 가중치를 곱해 최종 toxicity (각 0~12) 산출"""
    w = abx_risk[abx]
    age_tox    = min(12, round(base["age"]    * w["age"]))
    renal_tox  = min(12, round(base["renal"]  * w["renal"]))
    hepatic_tox= min(12, round(base["hepatic"]* w["hepatic"]))
    return {"age": age_tox, "renal": renal_tox, "hepatic": hepatic_tox}

# === (기존) 상태 필터용 간단 스코어 ===
def get_status_score(patient):
    age = patient.get('age', 0)
    creat = patient.get('renal_function', {}).get('creatinine', 0)
    if   age <= 10:        age_score = 1
    elif age <= 30:        age_score = 2
    elif age <= 50:        age_score = 3
    elif age <= 60:        age_score = 4
    else:                  age_score = 5
    if   creat <= 0.9:     creat_score = 1
    elif creat <= 1.2:     creat_score = 2
    elif creat <= 1.8:     creat_score = 3
    elif creat <= 2.3:     creat_score = 4
    else:                  creat_score = 5
    return age_score, creat_score

gram_nodes = ['gram_positive', 'gram_negative']
state_nodes = ['extreme_age', 'extreme_creatinine']
KG = nx.DiGraph()
KG.add_nodes_from(gram_nodes, bipartite=0)
KG.add_nodes_from(state_nodes, bipartite=2)
KG.add_nodes_from(abx_nodes, bipartite=1)

for abx, grams in abx_to_gram.items():
    for gram in grams:
        KG.add_edge(abx, gram, relation='effective_against')

# 상태(극단값)용 간단 엣지
for abx, risk in abx_risk.items():
    if risk['age'] >= 1.2:         # 가중치로 간접 판정
        KG.add_edge(abx, 'extreme_age', relation='is_toxic_to')
    if risk['renal'] >= 1.2:
        KG.add_edge(abx, 'extreme_creatinine', relation='is_toxic_to')

def get_patient_states(patient):
    age_score, creat_score = get_status_score(patient)
    states = []
    if age_score >= 4:
        states.append("extreme_age")
    if creat_score >= 4:
        states.append("extreme_creatinine")
    return states

def recommend_antibiotics(patient):
    gram = patient['gram_status']
    agent = patient['infectious_agent']
    allergy = set(patient.get('allergy', []))
    drug_inter = set(patient.get('drug_interactions', []))
    suscept = patient.get('susceptibility', {}).get(agent, {})
    log = []

    # 0. 환자 상태 요약
    patient_states = get_patient_states(patient)
    state_str = ", ".join(patient_states) if patient_states else "없음"
    log.append("━━━━━━━━━━━━━━━━━━━━━━━")
    log.append(f"🔍 극도의 신기능 저하/고령 필터링: {state_str}")
    log.append("━━━━━━━━━━━━━━━━━━━━━━━\n")

    # 1단계: Gram + 상태 기반 후보
    candidates = []
    eliminated_state = []
    for abx in abx_nodes:
        if not KG.has_edge(abx, gram):
            continue
        toxic = False
        toxic_reasons = []
        for state in patient_states:
            if KG.has_edge(abx, state) and KG[abx][state].get('relation') == 'is_toxic_to':
                toxic = True
                toxic_reasons.append(state)
        if toxic:
            eliminated_state.append(f"  · {abx} (is_toxic_to: {', '.join(toxic_reasons)})")
        else:
            candidates.append(abx)
    log.append("🔹 1단계: Gram+상태 기반 후보\n" + ("    " + ", ".join(candidates) if candidates else "    없음"))
    if eliminated_state:
        log.append("  ⮩ [제외 항목]\n" + "\n".join(eliminated_state))
    log.append("")

    # 2단계: 알러지/상호작용
    filtered = []
    eliminated = []
    reason_dict = {}
    for abx in candidates:
        reasons = []
        if abx in allergy:
            reasons.append("알러지")
        if abx in drug_inter:
            reasons.append("Drug Interaction")
        if reasons:
            eliminated.append(abx)
            reason_dict[abx] = ', '.join(reasons)
        else:
            filtered.append(abx)
    log.append("🔹 2단계: 알러지/Drug Interaction")
    if eliminated:
        for abx in eliminated:
            log.append(f"  · {abx}: {reason_dict[abx]}")
    else:
        log.append("    제외 없음")
    log.append("")

    # 3단계: 감수성(S)만 남김
    filtered2 = []
    suspi_excluded = []
    for abx in filtered:
        if suscept.get(abx, None) == "S":
            filtered2.append(abx)
        else:
            suspi_excluded.append(abx)
    log.append("🔹 3단계: 감수성(S)만 통과")
    if suspi_excluded:
        for abx in suspi_excluded:
            log.append(f"  · {abx}: 감수성 미달 (S 아님)")
    else:
        log.append("    제외 없음")
    log.append("")

    # 4단계: 최종 후보 나열
    log.append("🔹 4단계: 최종 후보")
    if filtered2:
        for abx in filtered2:
            log.append(f"  · {abx}")
    else:
        log.append("    (추천 항생제 없음)")
    log.append("")

    # ============== Toxicity Score (0~12) & TOPSIS ==============
    if filtered2:
        base = patient_baseline_toxicity(patient)

        # 항목별 toxicity 테이블 로그
        tox_info = ["━━━━━━━━━━━━━━━━━━━━━━━", "💊 [항생제별 Toxicity Score (0~12)]  (낮을수록 좋음)"]
        abx_tox_map = {}
        for abx in filtered2:
            tox = abx_specific_toxicity(abx, base)
            abx_tox_map[abx] = tox
            tox_info.append(
                f"  · {abx:12}: Age={tox['age']:>2} / Renal={tox['renal']:>2} / Hepatic={tox['hepatic']:>2}"
            )
        log += tox_info + [""]

        # TOPSIS 계산 (이상해 = [0,0,0], 최악해 = [12,12,12])
        data = np.array([[abx_tox_map[a]['age'], abx_tox_map[a]['renal'], abx_tox_map[a]['hepatic']] for a in filtered2])
        ideal = np.array([0.0, 0.0, 0.0])
        worst = np.array([12.0, 12.0, 12.0])

        d_plus  = np.linalg.norm(data - ideal, axis=1)   # 이상해와의 거리 (가까울수록 좋음)
        d_minus = np.linalg.norm(data - worst, axis=1)   # 최악해와의 거리 (멀수록 좋음)
        Ci = d_minus / (d_plus + d_minus + 1e-9)

        sorted_idx = np.argsort(-Ci)
        topsis_result = [f"{filtered2[i]} (Ci={Ci[i]:.3f})" for i in sorted_idx]

        log.append("⭐ [TOPSIS 기반 순위추천] (이상해=0,0,0 / 최악해=12,12,12)")
        for rec in topsis_result:
            log.append(f"  · {rec}")
    else:
        log.append("⭐ 추천항생제 없음 알아서 결정")

    log.append("━━━━━━━━━━━━━━━━━━━━━━━")
    return filtered2, log

def draw_kg():
    plt.figure(figsize=(10, 4))
    pos = {}
    for i, gram in enumerate(gram_nodes):
        pos[gram] = (0, i)
    for i, abx in enumerate(abx_nodes):
        pos[abx] = (1, i)
    for i, state in enumerate(state_nodes):
        pos[state] = (2, i*0.7 - 1)
    for node in KG.nodes():
        if node not in pos:
            pos[node] = (1, len(pos))
    node_colors = []
    for n in KG.nodes():
        if n in gram_nodes:
            node_colors.append("skyblue")
        elif n in state_nodes:
            node_colors.append("violet")
        else:
            node_colors.append("orange")
    nx.draw(KG, pos, with_labels=True, arrows=True, node_color=node_colors, node_size=1800, font_size=10)
    edge_labels = nx.get_edge_attributes(KG, 'relation')
    nx.draw_networkx_edge_labels(KG, pos, edge_labels=edge_labels, font_color="gray", font_size=9)
    plt.title("Antibiotic ↔ Gram/State Knowledge Graph")
    plt.axis('off')
    plt.tight_layout()
    return plt

# ---- Streamlit 인터페이스 ----
st.title("↓ 기본환자정보 및 항생제 감수성 정보를 통한 항생제 추천 (TOPSIS+Toxicity 0~12) ↓")

# 환자 선택
patient_idx = st.selectbox(
    "환자 선택",
    range(len(patients)),
    format_func=lambda i: f"{patients[i]['patient_id']} ({patients[i]['infectious_agent']})"
)
patient = patients[patient_idx]

# 환자 정보 표시 - 표 + 감수성 테이블
summary = {
    "ID": patient['patient_id'],
    "연령": patient['age'],
    "신장수치(Cr)": patient['renal_function'].get('creatinine', None),
    "빌리루빈": patient['hepatic_function'].get('bilirubin', None),
    "감염중증도": patient['infection_severity'],
    "감염균": patient['infectious_agent'],
    "Gram": patient['gram_status'],
    "중성구감소증": '있음' if patient.get('neutropenia') else '없음',
    "알러지": ", ".join(patient.get('allergy', [])) if patient.get('allergy') else "-",
    "Drug Interaction": ", ".join(patient.get('drug_interactions', [])) if patient.get('drug_interactions') else "-",
}
st.subheader("환자 정보 요약")
st.table(pd.DataFrame([summary]))

# 항생제별 감수성 표
st.subheader("항생제별 감수성")
abx_sir = patient['susceptibility'][patient['infectious_agent']]
df_sir = pd.DataFrame(list(abx_sir.items()), columns=["항생제", "SIR"])
st.dataframe(df_sir)

# 추천 결과/Reasoning Log
if st.button("항생제 추천/결과 보기"):
    result, log = recommend_antibiotics(patient)
    st.subheader("추천 항생제 (후보)")
    if result:
        for abx in result:
            st.markdown(f"- 💊 **{abx}**")
    else:
        st.warning("추천 항생제가 없습니다.")

    st.subheader("추천 Reasoning Log")
    st.text("\n".join(log))
