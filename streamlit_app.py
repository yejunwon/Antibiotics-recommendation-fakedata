import streamlit as st
import json
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import json
import os
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
    "cefTRIAXone sod 2g": ["gram_negative"],                   # 일부 GP 커버 가능하나 주 대상은 GN
    "cefePIMe 1g": ["gram_negative"],                          # GN 위주, P.aeruginosa 포함
    "Amoxclan duo(R) 437.5mg/62.5mg": ["gram_positive", "gram_negative"],  # 광범위, E.faecalis 포함
    "Meropenem 500mg": ["gram_positive", "gram_negative"]      # 광범위, ESBL 포함
}

# 나이·신기능·간기능 관련 위험 점수 (임상적 감안, 0~5 범위)
abx_risk = {
    "Tazoferan(R) 4.5g": {"age": 3, "creatinine": 4, "hepatic": 2},
    "cefaZOLin 1g": {"age": 2, "creatinine": 3, "hepatic": 1},
    "Azithromycin 250mg": {"age": 1, "creatinine": 0, "hepatic": 3},
    "cefTRIAXone sod 2g": {"age": 2, "creatinine": 1, "hepatic": 3},
    "cefePIMe 1g": {"age": 3, "creatinine": 5, "hepatic": 1},
    "Amoxclan duo(R) 437.5mg/62.5mg": {"age": 1, "creatinine": 3, "hepatic": 4},
    "Meropenem 500mg": {"age": 2, "creatinine": 4, "hepatic": 2}
}

def get_status_score(patient):
    age = patient.get('age', 0)
    creat = patient.get('renal_function', {}).get('creatinine', 0)
    bili = patient.get('hepatic_function', {}).get('bilirubin', 0)
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
    if   bili <= 0.9:      hepatic_score = 1
    elif bili <= 1.2:      hepatic_score = 2
    elif bili <= 1.8:      hepatic_score = 3
    elif bili <= 2.3:      hepatic_score = 4
    else:                  hepatic_score = 5
    return age_score, creat_score, hepatic_score

gram_nodes = ['gram_positive', 'gram_negative']
state_nodes = ['extreme_age', 'extreme_creatinine', 'extreme_hepatic']
KG = nx.DiGraph()
KG.add_nodes_from(gram_nodes, bipartite=0)
KG.add_nodes_from(state_nodes, bipartite=2)
KG.add_nodes_from(abx_nodes, bipartite=1)
for abx, grams in abx_to_gram.items():
    for gram in grams:
        KG.add_edge(abx, gram, relation='effective_against')
for abx, risk in abx_risk.items():
    if risk['age'] >= 4:
        KG.add_edge(abx, 'extreme_age', relation='is_toxic_to')
    if risk['creatinine'] >= 4:
        KG.add_edge(abx, 'extreme_creatinine', relation='is_toxic_to')
    if risk['hepatic'] >= 4:
        KG.add_edge(abx, 'extreme_hepatic', relation='is_toxic_to')

def get_patient_states(patient):
    age_score, creat_score, hepatic_score = get_status_score(patient)
    states = []
    if age_score >= 4:
        states.append("extreme_age")
    if creat_score >= 4:
        states.append("extreme_creatinine")
    if hepatic_score >= 4:
        states.append("extreme_hepatic")
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
    log.append(f"🔍 극도의 신기능 저하/고령/간기능 저하 필터링: {state_str}")
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

    # 4단계: 최종 추천
    log.append("🔹 4단계: 최종 추천")
    if filtered2:
        for abx in filtered2:
            log.append(f"  · {abx}")
    else:
        log.append("    (추천 항생제 없음)")
    log.append("")

    # Toxicity Score 요약 (항생제별, scaled to 0~12)
    age_score, creat_score, hepatic_score = get_status_score(patient)
    tox_info = ["━━━━━━━━━━━━━━━━━━━━━━━", "💊 [항생제 Toxicity Score]"]
    for abx in filtered2:
        a_risk = abx_risk[abx]['age']
        c_risk = abx_risk[abx]['creatinine']
        h_risk = abx_risk[abx]['hepatic']
        age_raw = age_score * a_risk
        creat_raw = creat_score * c_risk
        hepatic_raw = hepatic_score * h_risk
        age_tox = age_raw * 12 / 25.0
        creat_tox = creat_raw * 12 / 25.0
        hepatic_tox = hepatic_raw * 12 / 25.0
        tox_info.append(f"  · {abx:12}: age({age_score}×{a_risk})={age_raw}→{age_tox:.1f}  /  cr({creat_score}×{c_risk})={creat_raw}→{creat_tox:.1f}  /  hep({hepatic_score}×{h_risk})={hepatic_raw}→{hepatic_tox:.1f}")
    log += tox_info + [""]

    if filtered2:
        age_score, creat_score, hepatic_score = get_status_score(patient)
        A_list = []
        C_list = []
        H_list = []
        for abx in filtered2:
            a_risk = abx_risk[abx]['age']
            c_risk = abx_risk[abx]['creatinine']
            h_risk = abx_risk[abx]['hepatic']
            A_list.append((age_score * a_risk) * 12 / 25.0)
            C_list.append((creat_score * c_risk) * 12 / 25.0)
            H_list.append((hepatic_score * h_risk) * 12 / 25.0)
        data = np.array(list(zip(A_list, C_list, H_list)))
        ideal = np.array([0, 0, 0])
        anti_ideal = np.array([12, 12, 12])
        dist_to_ideal = np.linalg.norm(data - ideal, axis=1)
        dist_to_anti = np.linalg.norm(data - anti_ideal, axis=1)
        Ci = dist_to_anti / (dist_to_ideal + dist_to_anti + 1e-9)
        sorted_idx = np.argsort(-Ci)
        topsis_result = [f"{filtered2[i]} (Ci={Ci[i]:.3f})" for i in sorted_idx]
        log.append("⭐ [사용가능 항생제 순위추천 (TOPSIS)]")
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
st.title("↓ 기본환자정보 및 항생제 감수성 정보를 통한 항생제 추천 ↓")

# 환자 선택
patient_idx = st.selectbox(
    "환자 선택",
    range(len(patients)),
    format_func=lambda i: f"{patients[i]['patient_id']} ({patients[i]['infectious_agent']})"
)
patient = patients[patient_idx]

# 환자 정보 표시 - 표 + 감수성 테이블
import pandas as pd

# 1. 환자 주요정보 표로 요약
summary = {
    "ID": patient['patient_id'],
    "연령": patient['age'],
    "신장수치": patient['renal_function']['creatinine'],
    "빌리루빈": patient['hepatic_function']['bilirubin'],
    "감염중증도": patient['infection_severity'],
    "감염균": patient['infectious_agent'],
    "Gram": patient['gram_status'],
    "중성구감소증": '있음' if patient['neutropenia'] else '없음',
    "알러지": ", ".join(patient['allergy']) if patient['allergy'] else "-",
    "Drug Interaction": ", ".join(patient['drug_interactions']) if patient['drug_interactions'] else "-",
}
st.subheader("환자 정보 요약")
st.table(pd.DataFrame([summary]))

# 2. 항생제별 감수성 표
st.subheader("항생제별 감수성")
abx_sir = patient['susceptibility'][patient['infectious_agent']]
df_sir = pd.DataFrame(list(abx_sir.items()), columns=["항생제", "SIR"])
st.dataframe(df_sir)


# 추천 결과/Reasoning Log
if st.button("항생제 추천/결과 보기"):
    result, log = recommend_antibiotics(patient)
    st.subheader("추천 항생제")
    if result:
        for abx in result:
            st.markdown(f"- 💊 **{abx}**")
    else:
        st.warning("추천 항생제가 없습니다.")


    st.subheader("추천 Reasoning Log")
    st.text("\n".join(log))
