from dash import html

from kineticanalysis.utils.texts import t


def issue_footer():
    return html.P([t("misc.github_issue_prefix"),
                   html.A(t("misc.github_issue_link_text"), href=t("misc.github_issue_url")),
                   t("misc.github_issue_suffix"),
                   ])
