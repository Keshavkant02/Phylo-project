from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "soil_16s_tree_reading_printout.pdf"


def main() -> None:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import inch
    from reportlab.platypus import (
        PageBreak,
        Paragraph,
        Preformatted,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            name="Small",
            parent=styles["BodyText"],
            fontName="Helvetica",
            fontSize=9.2,
            leading=12,
            spaceAfter=4,
        )
    )
    styles.add(
        ParagraphStyle(
            name="Question",
            parent=styles["BodyText"],
            fontName="Helvetica-Bold",
            fontSize=10,
            leading=13,
            spaceBefore=6,
            spaceAfter=3,
        )
    )
    styles.add(
        ParagraphStyle(
            name="AnswerLine",
            parent=styles["BodyText"],
            fontName="Helvetica",
            fontSize=9.5,
            leading=14,
            textColor=colors.HexColor("#444444"),
        )
    )

    doc = SimpleDocTemplate(
        str(OUTPUT),
        pagesize=letter,
        rightMargin=0.62 * inch,
        leftMargin=0.62 * inch,
        topMargin=0.55 * inch,
        bottomMargin=0.55 * inch,
        title="Soil 16S Tree-Reading Printout",
    )

    story = []
    story.append(Paragraph("Soil 16S Tree-Reading Printout", styles["Title"]))
    story.append(Paragraph("Name: ____________________________    Date: ____________________________", styles["Small"]))
    story.append(Spacer(1, 0.08 * inch))
    story.append(
        Paragraph(
            "Big idea: a phylogenetic tree is read by ancestry. To decide which tips are most closely related, find the most recent common ancestor, not the tips that look closest on the page.",
            styles["Small"],
        )
    )
    story.append(Paragraph("Tips in this activity", styles["Heading2"]))
    tips = [
        ["Soil_ASV_A", "Bacillus_subtilis_168"],
        ["Soil_ASV_B", "Rhizobium_leguminosarum_IAM12609"],
        ["Pseudomonas_fluorescens_CCM2115", "Streptomyces_coelicolor_rrnD"],
        ["Acidobacterium_capsulatum_ATCC51196", ""],
    ]
    table = Table(tips, colWidths=[3.35 * inch, 3.35 * inch])
    table.setStyle(
        TableStyle(
            [
                ("FONT", (0, 0), (-1, -1), "Courier", 8.2),
                ("TEXTCOLOR", (0, 0), (-1, -1), colors.HexColor("#222222")),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#dddddd")),
                ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#fafafa")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ]
        )
    )
    story.append(table)
    story.append(Spacer(1, 0.1 * inch))

    story.append(Paragraph("Simplified Teaching Tree", styles["Heading2"]))
    tree = """root
|
+-- Acidobacterium_capsulatum_ATCC51196
|
+-- Streptomyces_coelicolor_rrnD
|
+-- Soil_ASV_A
|   |
|   +-- Bacillus_subtilis_168
|
+-- Pseudomonas_fluorescens_CCM2115
|
+-- Soil_ASV_B
    |
    +-- Rhizobium_leguminosarum_IAM12609"""
    story.append(Preformatted(tree, styles["Code"]))
    story.append(
        Paragraph(
            "For the full tree, use the notebook plots. This simplified drawing is for relationship-reading practice.",
            styles["Small"],
        )
    )
    story.append(Paragraph("Vocabulary", styles["Heading2"]))
    vocab = [
        ["Tip", "A sequence or organism shown at the end of a branch."],
        ["Node", "A branching point."],
        ["MRCA", "Most recent common ancestor."],
        ["Sister taxa", "Two tips or groups that share an immediate common ancestor."],
        ["Clade", "An ancestor and all of its descendants."],
        ["Branch length", "Sequence-distance information when drawn as a phylogram."],
    ]
    vocab_table = Table(vocab, colWidths=[1.25 * inch, 5.45 * inch])
    vocab_table.setStyle(
        TableStyle(
            [
                ("FONT", (0, 0), (0, -1), "Helvetica-Bold", 8.8),
                ("FONT", (1, 0), (1, -1), "Helvetica", 8.8),
                ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#dddddd")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 5),
                ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    story.append(vocab_table)

    story.append(PageBreak())
    story.append(Paragraph("Questions", styles["Heading1"]))
    questions = [
        "Which cached reference is Soil_ASV_A closest to?",
        "Which cached reference is Soil_ASV_B closest to?",
        "Is Bacillus_subtilis_168 the ancestor of Soil_ASV_A, or are both tips descendants of a shared ancestor?",
        "Which idea is used to decide relatedness on a tree? A. page closeness  B. color  C. most recent common ancestor  D. familiarity",
        "If a branch is rotated around a node, does the evolutionary relationship change? Explain.",
        "In this teaching tree, name one sister pair.",
        "Does a high 16S percent identity prove exact species identity? Explain.",
        "What does branch length show in this notebook's tree? A. p-value  B. sequence distance  C. abundance  D. moisture",
        "What statistic is more appropriate for support on a tree branch? A. adjusted p-value  B. bootstrap support  C. read count  D. sample depth",
        "Write one careful claim from this tree using: Soil_ASV_ ___ is closest to __________ in this cached 16S marker comparison, but this does not prove exact species identity.",
    ]
    for i, question in enumerate(questions, start=1):
        story.append(Paragraph(f"{i}. {question}", styles["Question"]))
        story.append(Paragraph("Answer: ________________________________________________________________", styles["AnswerLine"]))
        if i in {3, 5, 7, 10}:
            story.append(Paragraph("________________________________________________________________________", styles["AnswerLine"]))

    story.append(Spacer(1, 0.1 * inch))
    story.append(Paragraph("Exit Ticket", styles["Heading2"]))
    story.append(
        Paragraph(
            "In one sentence, explain the difference between identifying a closest reference and showing evolutionary relatedness.",
            styles["Small"],
        )
    )
    story.append(Paragraph("________________________________________________________________________", styles["AnswerLine"]))
    story.append(Paragraph("________________________________________________________________________", styles["AnswerLine"]))

    story.append(PageBreak())
    story.append(Paragraph("Instructor Answer Key", styles["Title"]))
    answers = [
        "Soil_ASV_A is closest to Bacillus_subtilis_168.",
        "Soil_ASV_B is closest to Rhizobium_leguminosarum_IAM12609.",
        "Both are tips descended from a shared ancestor. The living reference is not the ancestor of the ASV.",
        "C. Relatedness is decided by the most recent common ancestor.",
        "No. Rotation around a node changes the drawing, not the relationship.",
        "Accept Soil_ASV_A with Bacillus_subtilis_168, or Soil_ASV_B with Rhizobium_leguminosarum_IAM12609.",
        "No. High 16S identity supports a closest-reference claim, but one short 16S region does not prove exact species identity.",
        "B. Branch length is sequence-distance-derived branch length.",
        "B. Bootstrap support is more appropriate for branch support. Adjusted p-values belong to abundance or metadata association tests.",
        "Accept cautious claims such as: Soil_ASV_A is closest to Bacillus_subtilis_168 in this cached 16S marker comparison, but this does not prove exact species identity.",
    ]
    for i, answer in enumerate(answers, start=1):
        story.append(Paragraph(f"{i}. {answer}", styles["Small"]))

    story.append(Spacer(1, 0.15 * inch))
    story.append(Paragraph("Instructor Notes", styles["Heading2"]))
    notes = [
        "This printout is original and uses the soil 16S cache built for the Colab.",
        "It is inspired by tree-thinking skill targets from Baum, Smith, and Donovan, but does not reproduce their quiz pages.",
        "Emphasize: do not read trees as ladders of progress; do not decide relatedness by left-right order; do not treat tips as ancestors; do not call branch length a p-value.",
        "Use this before the Colab, then ask students to compare the simplified tree with the notebook's UPGMA and neighbor-joining plots.",
    ]
    for note in notes:
        story.append(Paragraph(f"- {note}", styles["Small"]))

    doc.build(story)
    print(OUTPUT)


if __name__ == "__main__":
    main()
