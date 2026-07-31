import streamlit as st
import traceback
import random
from src.template_loader import (
    get_branches, 
    get_areas, 
    get_template_files, 
    load_template_functions, 
    get_source_code
)
from src.storage import save_review

#  Page Configuration 
st.set_page_config(layout="wide", page_title="EngChain Annotator")

#  Custom CSS for Styling 
# This forces the primary buttons to have a specific look (optional but helps visibility)
st.markdown("""
<style>
    div.stButton > button:first-child {
        font-weight: bold;
    }
    /* FIX: Lock sidebar position and disable scrolling */
    section[data-testid="stSidebar"] {
        position: fixed;
        overflow: hidden !important;
    }
</style>
""", unsafe_allow_html=True)

#  Helper Function: Build the Queue 
def build_review_queue(branch):
    """
    Scans the entire branch and creates a list of all templates to review.
    Returns: List of dicts {'branch', 'area', 'file', 'name', 'func'}
    """
    queue = []
    areas = get_areas(branch)
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    total_steps = len(areas)
    for i, area in enumerate(areas):
        status_text.text(f"Loading area: {area}...")
        files = get_template_files(branch, area)
        for file in files:
            # Load all functions from this file
            try:
                funcs = load_template_functions(branch, area, file)
                for func_name, func_obj in funcs:
                    queue.append({
                        "branch": branch,
                        "area": area,
                        "file": file,
                        "name": func_name,
                        "func": func_obj
                    })
            except Exception as e:
                print(f"Error loading {file}: {e}")
        progress_bar.progress((i + 1) / total_steps)
        
    status_text.empty()
    progress_bar.empty()
    return queue

#  Session State Initialization 
if "app_mode" not in st.session_state:
    st.session_state["app_mode"] = "landing" # landing, review, done
if "review_queue" not in st.session_state:
    st.session_state["review_queue"] = []
if "current_index" not in st.session_state:
    st.session_state["current_index"] = 0
if "annotator_name" not in st.session_state:
    st.session_state["annotator_name"] = ""
if "current_branch" not in st.session_state:
    st.session_state["current_branch"] = ""
if "review_submitted" not in st.session_state:
    st.session_state["review_submitted"] = False

# ==========================================
# VIEW 1: LANDING PAGE
# ==========================================
if st.session_state["app_mode"] == "landing":
    st.title("🛡️ EngChain Verification Portal")
    
    st.markdown("""
    ### Welcome, Expert Annotator!
    Thank you for contributing to **EngChain**. Your task is to verify the correctness and quality 
    of our symbolic engineering templates.
    
    **Instructions:**
    1. Select your Engineering Domain.
    2. Enter your Name/ID.
    3. You will be guided through templates one-by-one.
    4. For each template, check the **Code** and the **Generated Trace**.
    5. Rate it, Approve/Reject it, and click **Submit**.
    """)
    
    st.markdown("")
    
    with st.form("onboarding_form"):
        st.subheader("Annotator Details")
        col1, col2 = st.columns(2)
        
        with col1:
            branches = get_branches()
            branch_input = st.selectbox("Select Your Domain", branches)
        
        with col2:
            name_input = st.text_input("Enter Your Name / ID")
            
        # Button styling: type="primary" gives it emphasis (color depends on theme, usually red/blue)
        submitted = st.form_submit_button("Start Annotation Session", type="primary")
        
        if submitted:
            if not name_input.strip():
                st.error("Please enter your name to proceed.")
            else:
                st.session_state["annotator_name"] = name_input
                st.session_state["current_branch"] = branch_input
                
                # Build the queue
                with st.spinner(f"Gathering templates for {branch_input}..."):
                    queue = build_review_queue(branch_input)
                
                st.session_state["review_queue"] = queue
                st.session_state["current_index"] = 0
                st.session_state["app_mode"] = "review"
                st.rerun()

# ==========================================
# VIEW 2: REVIEW WORKSPACE
# ==========================================
elif st.session_state["app_mode"] == "review":
    
    queue = st.session_state["review_queue"]
    idx = st.session_state["current_index"]
    
    # Check if we are done
    if idx >= len(queue):
        st.session_state["app_mode"] = "done"
        st.rerun()
    
    current_item = queue[idx]
    
    #  Sidebar Info 
    st.sidebar.title("Progress")
    st.sidebar.progress((idx) / len(queue))
    st.sidebar.write(f"Template: {idx + 1} / {len(queue)}")
    
    st.sidebar.markdown("")
    st.sidebar.subheader("Current Context")
    
    # Distinct styling for keys and values
    st.sidebar.markdown("**📂 Area:**")
    st.sidebar.info(current_item['area'])
    
    st.sidebar.markdown("**📄 File:**")
    st.sidebar.info(current_item['file'])
    
    st.sidebar.markdown("**🧩 Function:**")
    st.sidebar.info(current_item['name'])
    
    #  Main Content 
    st.title(f"Reviewing: {current_item['name']}")
    
    # Tabs for Inspection
    tab1, tab2 = st.tabs(["Source Code", "Generated Trace"])
    
    with tab1:
        code_text = get_source_code(current_item['func'])
        st.code(code_text, language="python", line_numbers=True)
        
    with tab2:
        col_gen, _ = st.columns([1, 4])
        # Renamed button for clarity
        if col_gen.button("Generate New Random Instance"):
            pass # Rerun trigger to get new random numbers
            
        try:
            question, solution = current_item['func']()
            st.markdown("#### Question")
            st.info(question)
            st.markdown("#### Solution Trace")
            st.success(solution)
        except Exception as e:
            st.error("Template Execution Failed")
            st.code(traceback.format_exc())

    st.markdown("")

    #  Scoring Form 
    st.subheader("Audit Decision")
    
    # Use a container so we can disable the form after submission
    with st.container():
        # Check if already submitted for this specific template
        is_submitted = st.session_state.get("review_submitted", False)
        
        if not is_submitted:
            # FORM STATE
            with st.form("audit_form"):
                c1, c2, c3 = st.columns(3)
                phys = c1.slider("Physical Plausibility", 1, 5, 5)
                math = c2.slider("Mathematical Correctness", 1, 5, 5)
                ped = c3.slider("Pedagogical Clarity", 1, 5, 5)
                
                decision = st.radio("Certification:", ["Approve", "Reject"], horizontal=True)
                feedback = st.text_area("Feedback (Required for Rejection)", placeholder="Explain any errors found...")
                
                # Removed emoji from button
                submit_review = st.form_submit_button("Submit Review")
                
                if submit_review:
                    if decision == "Reject" and not feedback.strip():
                        st.error("Feedback is required for Rejection.")
                    else:
                        # Save to disk
                        save_review(
                            st.session_state["annotator_name"],
                            current_item["branch"],
                            current_item["area"],
                            current_item["name"],
                            [phys, math, ped],
                            decision,
                            feedback
                        )
                        st.success("Review Saved!")
                        st.session_state["review_submitted"] = True
                        st.rerun()
        else:
            # POST-SUBMISSION STATE (Show Next Button)
            st.success("Review recorded for this template.")
            
            # Removed emoji, kept primary type for visibility
            if st.button("Proceed to Next Template", type="primary"):
                st.session_state["current_index"] += 1
                st.session_state["review_submitted"] = False
                st.rerun()

# ==========================================
# VIEW 3: COMPLETION PAGE
# ==========================================
elif st.session_state["app_mode"] == "done":
    # Removed balloons()
    st.title("Session Complete")
    
    st.success(f"""
    **Thank you, {st.session_state['annotator_name']}!**
    
    You have successfully reviewed all **{len(st.session_state['review_queue'])}** templates 
    in the **{st.session_state['current_branch']}** domain.
    """)
    
    st.info("Your reviews have been saved securely. You may now close this tab.")
    
    if st.button("Start New Session"):
        st.session_state.clear()
        st.rerun()
