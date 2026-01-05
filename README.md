# EngTrace: A Symbolic Benchmark for Verifiable Multi-Step Reasoning in Engineering

## Overview

EngTrace is a symbolic benchmark designed to evaluate the multi-step reasoning capabilities of Large Language Models (LLMs) within specialized engineering domains. Unlike static benchmarks that rely on factual recall, EngTrace requires models to synthesize scientific principles, mathematical modeling, and physical constraints to solve complex problems.

## Benchmark Structure

The repository contains 90 symbolic templates systematically organized into three primary engineering branches, nine core domains, and twenty pedagogical areas:

* **Chemical Engineering:** Reaction Kinetics, Thermodynamics, and Transport Phenomena.
* **Electrical Engineering:** Digital Communications, Electromagnetics, and Signals & Systems.
* **Mechanical Engineering:** Fluid Mechanics, Mechanics of Materials, and Vibrations & Acoustics.

Through domain-aware parameterization, these templates generate 1,350 unique, contamination-resistant test cases that ensure models cannot rely on memorized data.

## Key Features

* **Physically Grounded:** Problems are derived from real-world material properties and authoritative data sources to ensure physical realism.
* **Verifiable Reasoning:** Each problem includes a rigorous, step-by-step reasoning chain, allowing for the verification of the entire solution process rather than just the final numerical answer.
* **AI Tribunal & Quality Assurance:** The benchmark utilizes a tiered verification protocol where frontier reasoning models and domain experts certify the admissibility and logical consistency of the templates.

## Modules

* **data/templates/:** The core symbolic templates for Chemical, Electrical, and Mechanical engineering.
* **ai_assisted_quality_assurance/:** Tools for the "AI Tribunal" process to filter and validate problem templates.
* **evaluation/:** Scripts for running model inference, parsing engineering-specific outputs, and calculating AI-human alignment.
* **calculation_scripts:** Utilities for computing Inter-Annotator Agreement (IAA) and statistical alignment between model outputs and expert solutions.

