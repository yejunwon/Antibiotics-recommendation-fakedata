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

abx_nodes = [
    "ceftriaxone", "cefepime", "piperacillin/tazobactam",
    "meropenem", "imipenem", "ertapenem",
    "gentamicin", "amikacin",
    "vancomycin",
    "ciprofloxacin", "levofloxacin"
]

abx_to_gram = {
    "ceftriaxone": ["gram_negative"],
    "cefepime": ["gram_negative"],
    "piperacillin/tazobactam": ["gram_positive", "gram_negative"],
    "meropenem": ["gram_positive", "gram_negative"],
    "imipenem": ["gram_positive", "gram_negative"],
    "ertapenem": ["gram_positive", "gram_negative"],
    "gentamicin": ["gram_positive", "gram_negative"],
    "amikacin": ["gram_positive", "gram_negative"],
    "vancomycin": ["gram_positive"],
    "ciprofloxacin": ["gram_positive", "gram_negative"],
    "levofloxacin": ["gram_positive", "gram_negative"]
}

abx_risk = {
    "ceftriaxone": {"age": 1, "creatinine": 1},
    "cefepime": {"age": 2, "creatinine": 3},
    "piperacillin/tazobactam": {"age": 2, "creatinine": 2},
    "meropenem": {"age": 1, "creatinine": 2},
    "imipenem": {"age": 2, "creatinine": 2},
    "ertapenem": {"age": 1, "creatinine": 2},
    "gentamicin": {"age": 3, "creatinine": 5},
    "amikacin": {"age": 5, "creatinine": 4},
    "vancomycin": {"age": 2, "creatinine": 4},
    "ciprofloxacin": {"age": 2, "creatinine": 1},
    "levofloxacin": {"age": 2, "creatinine": 1}
}

score2onehot = {
    1:1, 2:2, 3:3, 4:4, 5:5, 6:6, 8:7, 9:8,
    10:9, 12:10, 15:11, 16:12, 20:13, 25:14
}
def get_onehot(score):
    return score2onehot.get(score, 0)
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
for abx, risk in abx_risk.items():
    if risk['age'] >= 4:
        KG.add_edge(abx, 'extreme_age', relation='is_toxic_to')
    if risk['creatinine'] >= 4:
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
    # 1. 기본 정보 추출
    # ... (age, creat, gram, agent 등) ...

    log = []

    # ===============================
    # 0. 데이터 결측 상황 분기 처리
    # ===============================
    if not agent and gram and (age is not None) and (creat is not None):
        # [Case 1] 감염균 정보만 없음 (Gram/나이/신기능 모두 있음)
        log.append("⏳ 감염균 정보 없음 → broad-spectrum antibiotics 우선 적용")
        candidates = [abx for abx in abx_nodes if set(abx_to_gram[abx]) == {'gram_positive', 'gram_negative'}]
        # → 이후 toxic, allergy filtering 등 기존 코드(공통 파이프라인)로 이어짐

    elif (age is None) or (creat is None):
        # [Case 2] 나이 또는 creatinine 없음
        log.append("⏳ 나이/신기능 정보 불충분 → broad-spectrum antibiotics 모두 표시, filtering 미적용")
        candidates = [abx for abx in abx_nodes if set(abx_to_gram[abx]) == {'gram_positive', 'gram_negative'}]
        # → allergy, drug interaction만 적용, 나머지는 skip 후 return

        # allergy/drug interaction filtering만 적용
        filtered = []
        for abx in candidates:
            if abx in allergy or abx in drug_inter:
                continue
            filtered.append(abx)
        log.append("🔹 알러지/Drug Interaction 제외 완료")
        log.append(f"추천: {', '.join(filtered) if filtered else '없음'}")
        return filtered, log

    elif (not agent) and (not gram):
        # [Case 3] 감염균/Gram/나이/신기능 정보 전부 없음
        log.append("⏳ 모든 정보 없음 → broad-spectrum antibiotics 전부 추천 (알러지 제외)")
        candidates = [abx for abx in abx_nodes if set(abx_to_gram[abx]) == {'gram_positive', 'gram_negative'}]
        # allergy만 제외 후 return
        filtered = []
        for abx in candidates:
            if abx in allergy:
                continue
            filtered.append(abx)
        log.append("🔹 알러지 제외 완료")
        log.append(f"추천: {', '.join(filtered) if filtered else '없음'}")
        return filtered, log

    else:
        # 기존 recommend_antibiotics 전체 코드 (Gram/나이/신기능/agent 전부 있는 경우)
        # 아래에서 candidates = [] 로부터 기존 파이프라인 그대로!
        # (즉, 기존 코드 그대로 사용)
        # ...
        # 마지막에 return filtered2, log
        pass


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

    # 4단계: 최종 추천
    log.append("🔹 4단계: 최종 추천")
    if filtered2:
        for abx in filtered2:
            log.append(f"  · {abx}")
    else:
        log.append("    (추천 항생제 없음)")
    log.append("")

    # Toxicity Score 요약 (항생제별)
    age_score, creat_score = get_status_score(patient)
    tox_info = ["━━━━━━━━━━━━━━━━━━━━━━━", "💊 [항생제 Toxicity Score]"]
    for abx in filtered2:
        a_risk = abx_risk[abx]['age']
        c_risk = abx_risk[abx]['creatinine']
        age_tox = get_onehot(age_score * a_risk)
        creat_tox = get_onehot(creat_score * c_risk)
        tox_info.append(f"  · {abx:12}: age({age_score}×{a_risk})={age_score*a_risk}->{age_tox}  /  cr({creat_score}×{c_risk})={creat_score*c_risk}->{creat_tox}")
    log += tox_info + [""]

    # TOPSIS 종합 순위
    if filtered2:
        age_score, creat_score = get_status_score(patient)
        A_list = []
        C_list = []
        for abx in filtered2:
            a_risk = abx_risk[abx]['age']
            c_risk = abx_risk[abx]['creatinine']
            A_list.append(get_onehot(age_score * a_risk))
            C_list.append(get_onehot(creat_score * c_risk))
        data = np.array(list(zip(A_list, C_list)))
        ideal = data.min(axis=0)
        anti_ideal = data.max(axis=0)
        dist_to_ideal = np.linalg.norm(data - ideal, axis=1)
        dist_to_anti = np.linalg.norm(data - anti_ideal, axis=1)
        Ci = dist_to_anti / (dist_to_ideal + dist_to_anti + 1e-9)
        sorted_idx = np.argsort(-Ci)
        topsis_result = [f"{filtered2[i]} (Ci={Ci[i]:.3f})" for i in sorted_idx]
        log.append("⭐ [사용가능 항생제 순위추천]")
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
st.title("항생제 추천 시스템 (Streamlit Demo)")

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

# Knowledge Graph 시각화 버튼
with st.expander("Knowledge Graph 시각화"):
    fig = draw_kg()
    st.pyplot(fig)
    plt.close()  # Streamlit에서 자원 누수 방지

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



