#!/bin/sh

# sync in generic modules

rsync -avz ../snippets/mtrand/mtrand.h prog/
rsync -avz ../snippets/vct/vct.h prog/
rsync -avz ../snippets/vct/nat.h prog/
rsync -avz ../snippets/wl/wl.h prog/
rsync -avz ../snippets/md/mdutil.h prog/
rsync -avz ../snippets/mtrand/mtrand.js web/js
rsync -avz ../snippets/vct/vct.js web/js
rsync -avz ../snippets/vct/mat.js web/js
rsync -avz ../snippets/wl/wl.js web/js
rsync -avz ../snippets/md/mdutil.js web/js
