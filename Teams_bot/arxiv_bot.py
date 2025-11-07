import feedparser, json, datetime, requests
import time
FLOW_URL = r"https://default723a5a87f39a4a2292473fc240c013.96.environment.api.powerplatform.com:443/powerautomate/automations/direct/workflows/90fa284934f1486f806f96b119ad7155/triggers/manual/paths/invoke?api-version=1&sp=%2Ftriggers%2Fmanual%2Frun&sv=1.0&sig=SUQbG59Bee34wkZ_gs9yReKiPp9NblwXmQYvDVR--jg" 

SEARCH_URL = (
    "http://export.arxiv.org/api/query?"
    "search_query=(ti:(qudit+AND+ion)+OR+abs:(qudit+AND+ion))"
    "&sortBy=submittedDate&sortOrder=descending&max_results=50"
)
SEEN_FILE = "seen.json"


def load_seen():
    try:
        with open(SEEN_FILE, "r") as f:
            return set(json.load(f))
    except:
        return set()


def save_seen(seen):
    with open(SEEN_FILE, "w") as f:
        json.dump(list(seen), f)


def fetch_latest():
    feed = feedparser.parse(SEARCH_URL)
    if not feed.entries:
        print("⚠️ No entries found.")
        return []
    return feed.entries


def check_new_papers():
    seen = load_seen()
    feed_entries = fetch_latest()
    new_papers = []

    for entry in feed_entries:
        if entry.id not in seen:
            seen.add(entry.id)
            new_papers.append(entry)

    save_seen(seen)
    return new_papers


def post_to_teams(title, link, summary, authors, date):
    payload = {
        "type": "message",
        "attachments": [
            {
                "contentType": "application/vnd.microsoft.card.adaptive",
                "content": {
                    "$schema": "http://adaptivecards.io/schemas/adaptive-card.json",
                    "type": "AdaptiveCard",
                    "version": "1.5",
                    "body": [
                        {
                            "type": "Container",
                            "items": [
                                {
                                    "type": "TextBlock",
                                    "text": title,
                                    "wrap": True,
                                    "weight": "Bolder",
                                    "size": "Large",
                                },
                                {
                                    "type": "TextBlock",
                                    "text": f"👨‍🔬 **Authors:** {authors}",
                                    "wrap": True,
                                    "spacing": "Small",
                                },
                                {
                                    "type": "TextBlock",
                                    "text": f"📅 **Published:** {date}",
                                    "wrap": True,
                                    "spacing": "Small",
                                },
                                {
                                    "type": "TextBlock",
                                    "text": "📝 **Abstract:**",
                                    "wrap": True,
                                    "weight": "Bolder"
                                },
                                {
                                    "type": "TextBlock",
                                    "text": summary,
                                    "wrap": True,
                                    "spacing": "Medium",
                                },
                            ],
                            "style": "emphasis",
                        }
                    ],
                    "actions": [
                        {
                            "type": "Action.OpenUrl",
                            "title": "🔗 Read on arXiv",
                            "url": link,
                        }
                    ],
                },
            }
        ],
    }

    r = requests.post(FLOW_URL, json=payload)
    print(f"→ Posted: {title} | HTTP {r.status_code}")


def summarize(entry):
    title = entry.title.strip()
    link = entry.link
    authors = ", ".join(a.name for a in entry.authors)
    date = entry.published.split("T")[0]
    summary = entry.summary.replace("\n", " ")
    return title, link, summary, authors, date


def main():
    print(f"\n🕓 Checking arXiv for new 'qudit trapped ion' papers (title+abstract)... ({datetime.date.today()})")
    new_papers = check_new_papers()

    if new_papers:
        print(f"✅ Found {len(new_papers)} new paper(s):\n")
        for e in new_papers[::-1]:
            title, link, summary, authors, date = summarize(e)
            print(f"- {title} | {authors} | {date}")
            post_to_teams(title, link, summary, authors, date)
            print("⏳ Waiting 1 minute before next post...")
            time.sleep(5)
    else:
        print("No new papers since last check.")


if __name__ == "__main__":
    main()
