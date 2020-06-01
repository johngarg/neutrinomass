#!/usr/bin/env python3

"""Functions to export networkx graphs to pdf with tikz-feynman."""

import os

from neutrinomass.completions.completions import operator_completions
from neutrinomass.completions.core import (
    FieldType,
    ComplexScalar,
    RealScalar,
    VectorLikeDiracFermion,
    MajoranaFermion,
    Completion,
)
from neutrinomass.completions.operators import EFF_OPERATORS
from neutrinomass.tensormethod.core import Operator, IndexedField

# g = operator_completions(EFF_OPERATORS["3b"])[0].graph

import networkx as nx


def skeleton(body: str) -> str:
    return (
        r"""
\documentclass{standalone}
\usepackage{tikz}
\usepackage[compat=1.1.0]{tikz-feynman}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{accents}

\begin{document}

   \feynmandiagram [spring electrical layout] {
       %s
   };

\end{document}
"""
        % body
    )


def tikz_type(particle: IndexedField) -> str:
    """Returns a string classifying the particle as (anti) fermion, (charged)
    scalar.

    """

    if isinstance(particle, ComplexScalar):
        return "charged scalar"

    if isinstance(particle, RealScalar):
        return "scalar"

    if isinstance(particle, MajoranaFermion) and not particle.is_conj:
        return "anti majorana"

    if isinstance(particle, MajoranaFermion) and particle.is_conj:
        return "majorana"

    if isinstance(particle, VectorLikeDiracFermion) and not particle.is_conj:
        return "anti majorana"

    if isinstance(particle, VectorLikeDiracFermion) and particle.is_conj:
        return "majorana"

    # SM particles
    if particle.is_scalar and particle.is_conj:  # Higgs
        return "anti charged scalar"

    if particle.is_scalar and not particle.is_conj:
        return "charged scalar"

    if particle.is_fermion and particle.is_conj:  # Quarks and leptons
        return "anti fermion"

    if particle.is_fermion and not particle.is_conj:
        return "fermion"

    return ValueError(f"Incomplete cases for {particle}.")


def proc_body(g: nx.Graph) -> str:
    external_nodes = [n for n in g.nodes if len(list(g.neighbors(n))) == 1]
    lines = []
    for a, b, data in g.edges.data():
        # edges that come about from
        # neutrinomass.completions.completion.get_connecting_edge are already
        # in the right order. The ones that aren't are the external edges,
        # these need to be ordered
        if a in external_nodes or b in external_nodes:
            a, b = sorted([a, b])

        p = data["particle"]
        fieldtype = tikz_type(p)
        lines.append(f"v{a} -- [{fieldtype}, edge label=${p.latex}$] v{b},")

    return "\n".join(lines)


def tikz_export(g: nx.Graph, path: str) -> None:
    """Exports latex document to path and builds output in texbuild directory.

    Example:
       >>> tikz_export(g, '~/Desktop/test.tex')

    """

    latex = skeleton(proc_body(g))

    # make directory
    file_dir, file_name = os.path.split(path)
    build_dir = os.path.join(file_dir, "texbuild")
    print(build_dir)
    mkdir = f"mkdir -p {build_dir}"
    os.system(mkdir)

    with open(path, "w+") as tex_file:
        tex_file.write(latex)

    # compile latex with lualatex
    compile_command = rf"lualatex -output-directory={build_dir} {path}; "
    os.system(compile_command)
    print(f"Successfully generated diagram at {path}!")


def tikz_export_generic(g: nx.Graph, path: str) -> None:
    pass


def export_completion(comp: Completion, path: str) -> None:
    pass
